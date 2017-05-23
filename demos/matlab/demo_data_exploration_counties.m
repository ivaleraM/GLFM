%% --------------------------------------------------
% DEMO: Data exploration on counties database
%% --------------------------------------------------
close all;
addpath(genpath('../../src/'));
rng(round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/counties.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias

data.X(:,9) = data.X(:,8) + data.X(:,9);
data.ylabel{9} = 'age >= 65';
% [1 3 4] remove dimensions with excessive number of missings
idx_to_remove = [1,3,4, 6, 7,8, 10,19]; % [1, 3, 4, 19]
data.X(:,idx_to_remove) = [];
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];

Xtrue = data.X;
% specify external transforms for certain dimensions
idx_transform = [2 5 6 11] ; %[4 5 10]; %[2 3 7 9 10 15];
params.t = cell(1,size(data.X,2));
params.t_1 = cell(1,size(data.X,2));
params.dt_1 = cell(1,size(data.X,2));
params.ext_dataType = cell(1,size(data.X,2));

for r=idx_transform
    params.t_1{r} = @(x) log(x + 1);
    params.t{r} = @(y) exp(y) - 1;
    params.dt_1{r} = @(x) 1 ./ (x + 1);
    params.ext_dataType{r} = 'p';
end
% % dimension 13 = 'White' need an inversion too
% r = 10; %9; %14;
% params.t_1{r} = @(x) log((100-x) + 1);
% params.t{r} = @(y) - exp(y) + 101;
% params.dt_1{r} = @(x) - 1./ (101 - x);
% params.ext_dataType{r} = 'p';

%% DEFINE PARAMS
[N, D] = size(data.X);

params.missing = -1;
%%params.s2Y = 0;     % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;  % Auxiliary variance
params.s2B = 1;     % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 100; % Number of iterations for the gibbs sampler
end
params.maxK = 10;

if ~isfield(params,'bias')
    params.bias = 1;
end
params.func = 1*ones(1,D);

if ~isfield(params,'simId')
    params.simId = 1;
end
if ~isfield(params,'save')
    params.save = 0;
end

saveFigs = 0;
savepath = './';

%% Initialize Hidden Structure

if (params.bias)
    Zini = [ones(N,1), double(rand(N,1)>0.8)];
else
    Zini = double(rand(N,1)>0.8);
end
hidden.Z = Zini; % N*D

%% Inference
hidden = GLFM_infer(data, hidden, params);

if params.save
    output_file = sprintf( './results/counties_bias%d_simId%d_Niter%d_s2B%.2f_alpha%d.mat', ...
        params.bias, params.simId, params.Niter, params.s2B, params.alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
if ~params.save
    X_map = GLFM_computeMAP(data.C, hidden.Z, hidden, params);
end

%% Plot Dimensions
if ~params.save

    sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*0.03);
    hidden = remove_dims(hidden, feat_toRemove);
    hidden = sort_hidden(hidden);
    sum(hidden.Z)

    [patterns, CC, L] = get_feature_patterns_sorted(hidden.Z);

    idxD = 2:size(data.X,2);

%     %Zp = patterns;
%     %leg = {'(1 0 0)', '(1 0 1)', '(1 1 0)', '(1 1 1)'};
%     Zp = eye(size(hidden.Z,2));
%     leg = num2str(Zp);
%     leg = mat2cell(leg, ones(size(hidden.Z,2),1), size(leg,2))';
    
    Zp = patterns(L > 240,:);
    leg = num2str(Zp(:,2:end));
    leg = mat2cell(leg, ones(size(Zp,1),1), size(leg,2))';
    %Zp = [Zp zeros(size(Zp,1), size(hidden.Z,2) - idxF)];
    
    colors = [ 0 102 255; ...
        153 51  255; ...
        204 204 0; ...
        255 102  102; ...
        0   204 102;
        255 51 255];
    colors = colors ./ 255; %repmat(sum(colors,2),1,3);
    colors(3,:) = [0.9290    0.6940    0.1250];
    colors(5,:) = [0.4660    0.6740    0.1880];
    
    % change order of colors
    colors = colors([3 5 4 2 1],:);
    colors(4,:) = [0 255 255] ./255;
    
    GLFM_plotPatterns(data, hidden, params, Zp, ...
        'leg', leg, 'colors', colors(1:size(Zp,1),:), 'idxD', idxD);
    if saveFigs
        cleanfigure;
        matlab2tikz([savepath, sprintf('dim%d.tex', d)] );
        saveas(gca,[savepath, sprintf('dim%d.fig', d)] );
        figurapdf(7,7)
        print([savepath, sprintf('dim%d.pdf', d)],['-d','pdf'],'-r150');
    end
end

%% PLOT USA map and corresponding features
if ~params.save

    for k=1:size(patterns,1)
        pat = patterns(k,:);
        Zn = (CC == k);
        %Zn = hidden.Z(:,k);

        if (sum(Zn) < (size(data.X,1)*0.01))
            continue; % only plot patterns with more than 5% of number of obs.
        end
        plot_usa_map(data,Zn);
        title(sprintf('Activation of pattern (%s)',num2str(pat)));
        if saveFigs
            %title(sprintf('Activation of pattern (%s)',num2str(pat)));
            cleanfigure;
            matlab2tikz([savepath, sprintf('clusters_map%d.tex', k)]);
            saveas(gca,[savepath, sprintf('clusters_map%d.fig', k) ] );
            figurapdf(7,7);
            print([savepath, sprintf('clusters_map%d.pdf', k) ],['-d','pdf'],'-r150');
        end
    end

end
