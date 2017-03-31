%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------
function gts_counties(s2B,func,simId,Niter, cluster)

if cluster
    output_folder = '/export/gts_usuarios/melanie/GLFM_results/';
else
    output_folder = './results/';
end

addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/counties.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias

% [1 3 4] remove dimensions with excessive number of missings
% dimension 14 = 'White' would need an inversion too, let's remove it for
% the moment
idx_to_remove = [1,3,4,7, 17];
data.X(:,idx_to_remove) = []; 
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];

Xtrue = data.X;
% pre-transform pop. variables
idx_transform = [2 3 6 9 13]; %7 9 10 15];
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

% %% NO PERC
% %
% mask_noMiss = ~isnan( data.X(:,14) );
% dem_id = (data.X(mask_noMiss,14) > data.X(mask_noMiss,15)) & (data.X(mask_noMiss,14) > data.X(mask_noMiss,16));
% rep_id = (data.X(mask_noMiss,15) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,15) > data.X(mask_noMiss,16));
% per_id = (data.X(mask_noMiss,16) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,16) > data.X(mask_noMiss,15));
% %
% idx_to_remove = [1,3,4, 14, 15, 16];
% data.X(:,idx_to_remove) = []; % remove dimensions with excessive number of missings
% data.C(idx_to_remove) = [];
% data.cat_labels(idx_to_remove) = [];
% data.ylabel(idx_to_remove) = [];

%% Initialize Hidden Structure

[N, D] = size(data.X);
%Zini = [dem_id, rep_id, per_id, double(rand(N,1)>0.8)];
Zini = [ones(N,1), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% DEFINE PARAMS
params.missing = -1;
params.s2Y = 0;     % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;  % Auxiliary variance
params.s2B = s2B;   % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = Niter; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
params.bias = 1;
params.func = func*ones(1,D);
params.simId = simId;
params.save = 1;

%params.simId = 1;
if ~isfield(params,'save')
    params.save = 1;
end

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = [ output_folder, sprintf( 'counties_bias%d_alpha%d_simId%d_Niter%d_s2B%.2f_func%d.mat', ...
        params.bias, params.alpha, params.simId, params.Niter, params.s2B, func) ];
    save(output_file);
end

%% Predict MAP estimate for each latent feature
if ~params.save
    X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);
end

%% PLOT USA map and corresponding features
if ~params.save
    
    for k=2:size(hidden.Z,2)
        
        plot_usa_map(data,hidden,k);
        title(sprintf('Activation of F%d',k-1));
        
        %     Zp = zeros(2,size(hidden.Z,2));
        %     Zp(:,1) = 1; % bias
        %     Zp(1,k) = 1; % feature active
        %
        %     X_F = IBPsampler_MAP(data.C, Zp, hidden);
        %     idx_toRemove = [1 4 14 15 16];
        %     X_F(:,idx_toRemove) = [];
        %     labels = data.ylabel;
        %     labels(idx_toRemove) = [];
        %
        %     plot_cont_feat(X_F, labels, sprintf('F%d',k-1));
    end
    
    plot_cont_all_feats(data, hidden, params);
end

%% Plot Dimensions
if ~params.save
    data.ylabel_long = data.ylabel;
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(6,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
    
    
    figure;
    for d=4:D
        subplot(2,1,1);
        [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            bar(pdf');
        elseif (data.C(d) == 'n')
            stem(xd, pdf');
        else
            plot(xd,pdf');
        end
        title(data.ylabel_long{d});
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            set(gca,'XTickLabel',data.cat_labels{d});
            set(gca,'XTickLabelRotation',45);
        end
        legend(leg);
        subplot(2,1,2);
        if ~isempty(params.t{d})
            data.X(:,d) = params.t_1{d}(data.X(:,d));
        end
        hist(data.X(:,d),100); title('Empirical');
        %     subplot(3,1,2);
        %     hist(data.X(drug_identifier,d),100); title('Empirical Bias 0');
        %     subplot(3,1,3);
        %     hist(data.X(drug_identifier,d),100); title('Empirical Bias 1');
        
        pause;
    end
end