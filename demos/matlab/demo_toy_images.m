%% demo TOY IMAGES

clear
addpath(genpath('../../src/GLFMmatlab/'));
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% GENERATIVE MODEL
N = 1000;
s2x = 1;
[data,gT] = generate_toy_images(N,s2x);
D = size(data.X,2);

%% INITIALIZATION + PARAMETER SETTINGS

Zini = double(rand(N,1)>0.8);
hidden.Z = Zini;

% define params
params.missing = -1;
params.s2Y = .5;      % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .01;
params.s2B = 1;       % Variance of the Gaussian prior of the weigting matrices B
params.alpha = 1;     % Concentration parameter of the IBP
params.Niter = 1000; % Number of iterations for the gibbs sampler
params.maxK = 10;
params.bias = 0;

%% INFERENCE
hidden = IBPsampler_run(data, hidden, params);

%% PLOT EACH LATENT FEATURE
Kest = size(hidden.B,2);
Zp = eye(Kest);
X_map = IBPsampler_MAP(Zp,hidden);
