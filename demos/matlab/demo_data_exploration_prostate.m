
%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------

addpath(genpath('../../src/'));
rng( round(sum(1e5*clock)) );
savePath = './results/';

%% LOAD DATA
input_file = '../../datasets/mat/prostate_v3.mat';
load(input_file);

% Reduce dataset (we only focus on a few dimensions, since we have little num. of obs.

idx_toKeep = [1 2 4 13 15];
data.X = data.X(:,idx_toKeep);
data.C = data.C(idx_toKeep);
data.cat_labels = data.cat_labels(idx_toKeep);
data.ylabel = data.ylabel(idx_toKeep);
data.ylabel_long = data.ylabel_long(idx_toKeep);


%% DEFINE PARAMS (optional parameters for GLFM model)
[N, D] = size(data.X);

% pre-transform a subset of variables
idx_transform = D; % we transform the last dimension
params.t = cell(1,D);
params.t_1 = cell(1,D);
params.dt_1 = cell(1,D);
params.ext_dataType = cell(1,D);
for r=idx_transform
    params.t_1{r} = @(x) log(x + 1); % transformation to apply to raw data
    params.t{r} = @(y) exp(y) - 1;   % inverse transform to recover raw data
    params.dt_1{r} = @(x) 1./(x+1);  % derivative of inverse transform
    params.ext_dataType{r} = 'p';    % change type of data due to transformation
end


params.missing = -1;
params.s2Y = 0;         % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;      % Auxiliary variance
if ~isfield(params,'s2B')
    params.s2B = 1;     % Variance of the Gaussian prior of the weigting matrices B
end
params.alpha = 1;       % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 100; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
if ~isfield(params,'bias')
    params.bias = 1;    % 1 = fix first feature to be active for all patients 
end
params.func = ones(1,D);

if ~isfield(params,'simId')
    params.simId = 1;  % simulation identifier (to be able to run the sample sim. multiple times 
end
if ~isfield(params,'save')
    params.save = 0;   % save .mat if active
end

%% Initialize Hidden Structure

if params.bias
    Zini = [ones(N,1), double(rand(N,1)>0.8)];
else
    Zini = double(rand(N,1)>0.8);
end
hidden.Z = Zini; % N*K matrix of feature assignments

%% Inference
hidden = IBPsampler_infer(data, hidden, params);

if params.save
    output_file = [savePath, sprintf('prostateRed_bias%d_simId%d_Niter%d_s2Y%.2f_s2B%.2f_alpha%.2f.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha)];
    save(output_file);
end

%% Predict MAP estimate for the whole matrix X
X_map = IBPsampler_computeMAP(data.C, hidden.Z, hidden, params);

%% Plot Dimensions
if ~params.save
    
    sum(hidden.Z)
    th = 0.03; % threshold to filter out latent features that are not significant
    feat_toRemove = find(sum(hidden.Z) < N*th); % filter features with insufficient number of obs. assigned
    hidden = remove_dims(hidden, feat_toRemove); % remove latent dimensions
    
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    % choose patterns corresponding to activation of each feature
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(5,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5'};
    
    colors = []; styles = [];
    plot_all_dimensions(data, hidden, params, Zp, leg, colors, styles);
end