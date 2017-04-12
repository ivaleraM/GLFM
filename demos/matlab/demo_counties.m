%% --------------------------------------------------
% DEMO: Data exploration on counties database
%% --------------------------------------------------
addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/counties.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias

% [1 3 4] remove dimensions with excessive number of missings
idx_to_remove = [1,3,4,19];
data.X(:,idx_to_remove) = [];
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];

Xtrue = data.X;
% specify external transforms for certain dimensions
idx_transform = [2 3 7 9 10 15]; %7 9 10 15];
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
% dimension 13 = 'White' need an inversion too
r = 14;
params.t_1{r} = @(x) log((100-x) + 1);
params.t{r} = @(y) - exp(y) + 101;
params.dt_1{r} = @(x) - 1./ (101 - x);
params.ext_dataType{r} = 'p';

%% DEFINE PARAMS
[N, D] = size(data.X);

params.missing = -1;
params.s2Y = 0;     % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;  % Auxiliary variance
params.s2B = 1;     % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 10; % Number of iterations for the gibbs sampler
end
params.maxK = 10;

if ~isfield(params,'bias')
    params.bias = 1;
end
params.func = 1*ones(1,D);

%params.simId = 1;
if ~isfield(params,'save')
    params.save = 0;
end

%% Initialize Hidden Structure

if (params.bias)
    Zini = [ones(N,1), double(rand(N,1)>0.8)];
else
    Zini = [double(rand(N,2)>0.8)];
end
hidden.Z = Zini; % N*D

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = sprintf( './results/counties_bias%d_simId%d_Niter%d_s2Y%.2f_s2B%.2f_alpha%d.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
if ~params.save
    X_map = IBPsampler_MAP(data.C, hidden.Z, hidden, params);
end

%% Plot Dimensions

if ~params.save
    
    sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*0.03);
    hidden = remove_dims(hidden, feat_toRemove);
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    params.th = 0.1;
    idxD = 4:size(data.X,2);
    Zp = patterns;
    leg = {'(1 0 0)', '(1 0 1)', '(1 1 0)', '(1 1 1)'};
    plot_all_dimensions(data, hidden, params, Zp, leg, idxD);
end

%% PLOT USA map and corresponding features
if ~params.save
    
    for k=1:size(patterns,1)
        pat = patterns(k,:);
        Zn = (C == k);%hidden.Z(:,idxF);
        plot_usa_map(data,Zn);
        title(sprintf('Activation of pattern (%s)',num2str(pat)));
    end
    
    plot_cont_all_feats(data, hidden, params);
end