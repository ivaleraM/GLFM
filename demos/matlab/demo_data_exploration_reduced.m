%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------
%clear
addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/prostate.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias
data.cat_labels{1} = {'3';'4'};

% change to ordinal variables

% data.C(3) = 'o'; % months to follow up
% tmp = num2str( unique(data.X(~isnan(data.X(:,3)),3)) );
% data.cat_labels{3} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));

% Type of activity: confined in bed, in bed less than 1/2 daytime, ...
data.C(7) = 'o';
tmp = (data.X(:,7) == 2);
data.X( data.X(:,7) == 3,7) = 2;
data.X(tmp,7) = 3;
tmp = data.cat_labels{7}(2);
data.cat_labels{7}(2) = data.cat_labels{7}(3);
data.cat_labels{7}(3) = tmp;
%tmp = data.cat_labels_long{7}(2);
%data.cat_labels_long{7}(2) = data.cat_labels_long{7}(3);
%ata.cat_labels_long{7}(3) = tmp;

% Index Stage (state of the patient)
data.C(14) = 'o';
tmp = num2str( unique(data.X(~isnan(data.X(:,14)),14)) );
data.cat_labels{14} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));

% ---------------------
% drug_identifier = data.X(:,2) > 0.5;
% % remove drug levels
% data.X(:,2) = [];
% data.C(2) = [];
% data.cat_labels(2) = [];
% data.ylabel(2) = [];
% data.ylabel_long(2) = [];

tmp = num2str(unique(data.X(:,2)));
data.cat_labels{2} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));
data.C(2) = 'o';

%% DEFINE PARAMS
[N, D] = size(data.X);

params.missing = -1;
params.s2Y = 0;       % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;      % Auxiliary variance
%params.s2B = 0.5;       % Variance of the Gaussian prior of the weigting matrices B
if ~isfield(params,'s2B')
    params.s2B = 1;   % Variance of the Gaussian prior of the weigting matrices B
end
params.alpha = 1;      % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 100; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
if ~isfield(params,'bias')
    params.bias = 1;
end
params.func = ones(1,D);

if ~isfield(params,'simId')
    params.simId = 1;
end
if ~isfield(params,'save')
    params.save = 0;
end

%% Reduced

% simplify Status
data.cat_labels{4};
% 1 -->      1  alive
% [2,3] -->  2  vascular
% 6 -->      3  prostatic cancer
% [7,8] -->  4  lung-related dead
% [4,5,9,10] 5  others
V = data.X(:,4);
V(data.X(:,4) == 3) = 2;
V(data.X(:,4) == 6) = 3;
%V(data.X(:,4) == 7 | data.X(:,4) == 8) = 4;
%V(data.X(:,4) == 4 | data.X(:,4) == 5 | data.X(:,4) == 9 | data.X(:,4) == 10) = 5;
V(data.X(:,4) == 7 | data.X(:,4) == 8 | ...
    data.X(:,4) == 4 | data.X(:,4) == 5 | data.X(:,4) == 9 | data.X(:,4) == 10) = 4;
%data.cat_labels{4} = {'alive', 'vascular', 'prostatic ca', 'lung-related', 'others'};
data.cat_labels{4} = {'alive', 'vascular', 'prostatic ca', 'others'};
data.X(:,4) = V;

% for drug level, fusion 0 and 0.2
mask = logical(data.X(:,2) < 0.5);
data.X(mask,2) = 0;
data.cat_labels{2}(2,:) = [];

idx_toKeep = [1 2 4 13 15]; % 12 ]; %[1 2 4 8 13 ]; %15]; %[1 2 3 4 5 8  13 15]; %[1 2 4 5 8 13 14];
%bias = data.X(:,1) - 3;
%hidden.Z = [bias, double(rand(N,1)>0.8)];
data.X = data.X(:,idx_toKeep);
data.C = data.C(idx_toKeep);
data.cat_labels = data.cat_labels(idx_toKeep);
data.ylabel = data.ylabel(idx_toKeep);
data.ylabel_long = data.ylabel_long(idx_toKeep);

drug_identifier = data.X(:,2) > 0.5;
drug_identifier_low = data.X(:,2) == 1; % high level drug
drug_identifier_high = data.X(:,2) == 5; % high level drug
%Zini = [ones(N,1), drug_identifier];
Zini = [ones(N,1), drug_identifier_high, drug_identifier_low]; %, double(rand(N,1)>0.9)];

% data.X(:,2) = [];
% data.C(2) = [];
% data.cat_labels(2) = [];
% data.ylabel(2) = [];
% data.ylabel_long(2) = [];

% pre-transform pop. variables
idx_transform = size(data.X,2); %7 9 10 15];
params.t = cell(1,size(data.X,2));
params.t_1 = cell(1,size(data.X,2));
params.dt_1 = cell(1,size(data.X,2));
params.ext_dataType = cell(1,size(data.X,2));
for r=idx_transform
    params.t_1{r} = @(x) log(x + 1);
    params.t{r} = @(y) exp(y) - 1;
    params.dt_1{r} = @(x) 1./(x+1);
    params.ext_dataType{r} = 'p';
end

%% Initialize Hidden Structure
%Zini = [drug_identifier, not(drug_identifier), double(rand(N,1)>0.8)];
Zini = [ones(N,1), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = sprintf( './results/prostateRed_bias%d_simId%d_Niter%d_s2Y%.2f_s2B%.2f_alpha%d.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
X_map = IBPsampler_MAP(data.C, hidden.Z, hidden, params);

%% Plot Dimensions
if ~params.save
    th = 0.03;
    %idxD = 1:size(data.X,2);
    
    sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*th);
    hidden = remove_dims(hidden, feat_toRemove);
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    % choose patterns corresponding to activation of each feature
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(5,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5'};
    
    plot_all_dimensions(data, hidden, params, Zp, leg);
end