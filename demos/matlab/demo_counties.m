%% --------------------------------------------------
% DEMO: Data exploration on counties database
%% --------------------------------------------------
clear
addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/counties.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias

% [1 3 4] remove dimensions with excessive number of missings
idx_to_remove = [1,3,4,7];
data.X(:,idx_to_remove) = [];
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];

Xtrue = data.X;
% pre-transform pop. variables
idx_transform = [2 3 6 9 14]; %7 9 10 15];
params.t = cell(1,size(data.X,2));
params.t_1 = cell(1,size(data.X,2));
params.dt_1 = cell(1,size(data.X,2));
for r=idx_transform
    params.t{r} = @(x) log(x + 1);
    params.t_1{r} = @(y) exp(y) - 1;
    params.dt_1{r} = @(y) exp(y);
    data.X(:,r) = params.t{r}(data.X(:,r)); % work in logarithm space better
    data.C(r) = 'p';
end
% dimension 13 = 'White' need an inversion too
r = 13;
params.t{r};
params.t{r} = @(x) log((100-x) + 1);
params.t_1{r} = @(y) - exp(y) + 101;
params.dt_1{r} = @(y) - exp(y);
data.X(:,r) = params.t{r}(data.X(:,r)); % work in logarithm space better
data.C(r) = 'p';

%% Initialize Hidden Structure

[N, D] = size(data.X);
Zini = [ones(N,1), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% DEFINE PARAMS
params.missing = -1;
params.s2Y = 1;     % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;  % Auxiliary variance
params.s2B = 1;     % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 1000; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
params.bias = 1;
params.func = 1*ones(1,D);

%params.simId = 1;
if ~isfield(params,'save')
    params.save = 0;
end

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = sprintf( './results/counties_bias%d_simId%d_Niter%d_s2Y%.2f_s2B%.2f_alpha%d.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
if ~params.save
    X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);
end

%% PLOT USA map and corresponding features
if ~params.save
    
    sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*0.03);
    hidden = remove_dims(hidden, feat_toRemove);
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    for k=1:size(patterns,1)
        pat = patterns(k,:);
        Zn = (C == k);%hidden.Z(:,idxF);
        plot_usa_map(data,Zn);
        title(sprintf('Activation of pattern (%s)',num2str(pat)));
    end
    
    plot_cont_all_feats(data, hidden, params);
end

data.ylabel_long = data.ylabel;

%% Plot Dimensions

    %         if ~isempty(params.t{d})
    %             data.X(:,d) = params.t_1{d}(data.X(:,d));
    %         end

if params.save
    params.th = 0.1;
    idxD = 4:size(data.X,2);
    Zp = patterns;
    leg = {'(1 0 0)', '(1 0 1)', '(1 1 0)', '(1 1 1)'};
    plot_all_dimensions(data, hidden, params, Zp, leg, idxD);
end
