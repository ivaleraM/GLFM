% -------------------
%% demo TOY IMAGES
% -------------------
clear
addpath( genpath('../../src/') );

rng( round(sum(1e5*clock)) );

%% GENERATIVE MODEL
N =1000;        % number of observations to generate
s2x = 0.5;      % observation noise
[data,gT] = generate_toy_images(N,s2x);

%% INITIALIZATION + PARAMETER SETTINGS
hidden.Z = double(rand(N,1)>0.8); % initialize N*K feature assignment matrix

% define params
params.missing = -10;%
% params.s2Y = 1;      % Variance of the Gaussian prior for the pseudo-observations
params.alpha = 5;    % Concentration parameter of the IBP
params.Niter = 100;  % Number of iterations for the Gibbs sampler
params.maxK = 10;    % Maximum number of latent features (for memory allocation)

%% INFERENCE
hidden = IBPsampler_infer(data, hidden, params);

%% PLOT REAL Vs INFERRED LATENT FEATURE
Kest = size(hidden.B,2);
Zp = eye(Kest);
% compute observations resulting when each of the latent features is active
X_map = IBPsampler_computeMAP(data.C, Zp, hidden, params);

figure(1); title('Ground truth');
for k=1:size(gT.B,1)
    subplot(3,3,k);
    imagesc(reshape(gT.B(k,:), [6 6]));
end
hold off;

figure(2); title('Inferred Latent Features');
for k=1:Kest
    subplot(3,3,k);
    imagesc(reshape(X_map(k,:), [6 6]));
end
hold off;
