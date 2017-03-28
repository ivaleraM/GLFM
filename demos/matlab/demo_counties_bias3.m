%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------
%clear
addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/counties.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias

% idx_to_remove = [1,3,4];
% data.X(:,idx_to_remove) = []; % remove dimensions with excessive number of missings
% data.C(idx_to_remove) = [];
% data.cat_labels(idx_to_remove) = [];
% data.ylabel(idx_to_remove) = [];

%% NO PERC
%
% for data exploration, remove nan
idx_unknownCounties = find( sum(isnan(data.X(:,14)),2) ); 
data.X(idx_unknownCounties,:) = [];

mask_noMiss = ~isnan( data.X(:,14) );
dem_id = (data.X(mask_noMiss,14) > data.X(mask_noMiss,15)) & (data.X(mask_noMiss,14) > data.X(mask_noMiss,16));
rep_id = (data.X(mask_noMiss,15) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,15) > data.X(mask_noMiss,16));
per_id = (data.X(mask_noMiss,16) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,16) > data.X(mask_noMiss,15));
%
idx_to_remove = [1,3,4, 14, 15, 16];
data.X(:,idx_to_remove) = []; % remove dimensions with excessive number of missings
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];

%% Initialize Hidden Structure

[N, D] = size(data.X);
Zini = [dem_id, rep_id, per_id, double(rand(N,1)>0.8)];
%Zini = [ones(N,1), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% DEFINE PARAMS
params.missing = -1;
params.s2Y = 1;     % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;  % Auxiliary variance
params.s2B = 0.5;   % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 1000; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
params.bias = 3;
params.func = 2*ones(1,D);

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
X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);

%% Plot Dimensions
if ~params.save
    data.ylabel_long = data.ylabel;

    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp = Zp(1:3,:);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    %Zp(:,1) = 1; % bias active
    %Zp = repmat(Zp,2,1); Zp(4:end,4) = 1;
    leg = {'F0','F1', 'F2', 'F0b',  'F1b',  'F2b'};
    
    
    figure(1);
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
        hist(data.X(:,d),100); title('Empirical');
        %     subplot(3,1,2);
        %     hist(data.X(drug_identifier,d),100); title('Empirical Bias 0');
        %     subplot(3,1,3);
        %     hist(data.X(drug_identifier,d),100); title('Empirical Bias 1');
        pause;
    end
end