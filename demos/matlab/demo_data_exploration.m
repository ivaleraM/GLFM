%% --------------------------------------------------
% DEMO: Data exploration on prostate cancer database
%% --------------------------------------------------

addpath(genpath('../../src/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD DATA
input_file = '../../datasets/mat/prostate.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias
drug_identifier = data.X(:,2) > 0.5;
% remove drug levels
data.X(:,2) = [];
data.C(2) = [];
data.cat_labels(2) = [];
data.ylabel(2) = [];
data.ylabel_long(2) = [];

%% Initialize Hidden Structure
[N, D] = size(data.X);
Zini = [drug_identifier, not(drug_identifier), double(rand(N,1)>0.8)];
hidden.Z = Zini; % N*D

%% DEFINE PARAMS
params.missing = -1;
params.s2Y = 1;   % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005;   % Auxiliary variance
params.s2B = 0.5;   % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;   % Concentration parameter of the IBP
params.Niter = 100;   % Number of iterations for the gibbs sampler
params.maxK = 10;
params.bias = 2;

params.simId = 1;
params.save = 0;

%% Inference
hidden = IBPsampler_run(data, hidden, params);

if params.save
    output_file = sprintf( '../results/prostate_simId%d_%d_s2Y%.2f_s2u%.2f_alpha%.2f.mat',simId,Niter, s2Y, s2u, alpha);
    save(output_file);
end

%% Predict MAP estimate for each latent feature
Kest = size(hidden.B,2);
Zp = eye(Kest);
Zp = Zp(1:2,:);
leg = {'Bias 0','Bias 1'};


X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);

%% Plot Dimensions
figure();
for d=1:D
    subplot(3,1,1);
    [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
    if (data.C(d) == 'c') 
        bar(pdf');
    elseif (data.C(d) == 'n') || (data.C(d) == 'o')
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
    subplot(3,1,2);
    hist(data.X(drug_identifier,d),100); title('Empirical Bias 0');
    subplot(3,1,3);
    hist(data.X(drug_identifier,d),100); title('Empirical Bias 1');
    pause;
end
