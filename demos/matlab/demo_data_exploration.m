% %% --------------------------------------------------
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

% Type of activity: confined in bed, in bed less than 1/2 daytime, ...
data.C(7) = 'o';
tmp = (data.X(:,7) == 2);
data.X( data.X(:,7) == 3,7) = 2;
data.X(tmp,7) = 3;
tmp = data.cat_labels{7}(2);
data.cat_labels{7}(2) = data.cat_labels{7}(3);
data.cat_labels{7}(3) = tmp;
%data.X( data.X(:,7) == 1), 7) == 2;

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
params.s2Y = 1;   % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;   % Auxiliary variance
<<<<<<< HEAD
params.s2B = 1;   % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
params.Niter = 100;   % Number of iterations for the gibbs sampler

=======
params.s2B = 0.5;   % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 10;   % Concentration parameter of the IBP
if ~isfield(params,'save')
    params.Niter = 100;   % Number of iterations for the gibbs sampler
end
>>>>>>> origin/master
params.maxK = 10;
params.bias = 1;
params.func = 2*ones(1,D);

params.simId = 1;
if ~isfield(params,'save')
     params.save = 0;
end

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = sprintf( './results/prostate_bias%d_simId%d_Niter%d_s2Y%.2f_s2B%.2f_alpha%d.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);

%% Plot Dimensions
if 1%params.save
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(3,Kest),:);
    leg = {'F0','F1', 'F2'};
    
    
    figure(1);
    for d=1:D
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
