%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------
function gts_prostate(s2B,func,simId,Niter, cluster)

if cluster
    output_folder = '/export/gts_usuarios/melanie/GLFM_results/';
else
    output_folder = './results/';
end

addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/prostate.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias
data.cat_labels{1} = {'3';'4'};

% change to ordinal variables

data.C(3) = 'o'; % months to follow up
tmp = num2str( unique(data.X(~isnan(data.X(:,3)),3)) );
data.cat_labels{3} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));

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


%% Initialize Hidden Structure
[N, D] = size(data.X);
%Zini = [drug_identifier, not(drug_identifier), double(rand(N,1)>0.8)];
Zini = [ones(N,1), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% DEFINE PARAMS
params.missing = -1;
params.s2Y = 0;       % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;      % Auxiliary variance
%params.s2B = 0.5;       % Variance of the Gaussian prior of the weigting matrices B
if ~isfield(params,'s2B')
params.s2B = s2B;   % Variance of the Gaussian prior of the weigting matrices B
end
params.alpha = 1;      % Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = Niter; % Number of iterations for the gibbs sampler
end
params.maxK = 10;
params.bias = 1;
params.func = func*ones(1,D);
params.simId = simId;

%params.simId = 1;
if ~isfield(params,'save')
    params.save = 1;
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

% idx_toKeep = [1 2 4 5 8 13 14]; % 11
% %bias = data.X(:,1) - 3;
% %hidden.Z = [bias, double(rand(N,1)>0.8)];
% data.X = data.X(:,idx_toKeep);
% data.C = data.C(idx_toKeep);
% data.cat_labels = data.cat_labels(idx_toKeep);
% data.ylabel = data.ylabel(idx_toKeep);
% data.ylabel_long = data.ylabel_long(idx_toKeep);

%% Inference
hidden = IBPsampler_run(data, hidden, params);
params.save = 1;

if params.save
    output_file = [ output_folder, sprintf( 'prostateRed_bias%d_alpha%d_simId%d_Niter%d_s2B%.2f_func%d.mat', ...
        params.bias, params.alpha, params.simId, params.Niter, params.s2B, func) ];
    save(output_file);
end

%% Predict MAP estimate for each latent feature
if ~params.save
    X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);
end

%% Plot Dimensions
if ~params.save
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(5,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5'};
    
    
    figure(1);
    for d=1:size(data.X,2)
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