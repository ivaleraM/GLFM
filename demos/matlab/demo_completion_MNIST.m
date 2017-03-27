%% demo TOY COMPLETE MNIST

clear
addpath(genpath('../../src/GLFMmatlab/'));
addpath(genpath('../../src/Ccode/'));

randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

%% LOAD MNIST dataset
% This DB can be downloaded here: http://yann.lecun.com/exdb/mnist/
images = loadMNISTImages('../../datasets/raw/MNIST/train-images-idx3-ubyte');
N = 50; % subset of images to consider
Xtrue = images(:, randperm(size(images,2),N))' + 1; % N*D matrix (each column is one pixel

perc_missing = 0.2;
missing_val = -100;
mask_missings = rand(size(Xtrue)) < perc_missing;
Xmiss = Xtrue;
Xmiss(mask_missings) = missing_val;

C = repmat('n',1,size(Xmiss,2));

data.X = Xmiss;
data.C = C;

%% Define parameter default values for algorithm
params.bias = 1;
params.s2Y = 1;
params.s2B = 1;
params.alpha = 10;
params.Niter = 100;
params.maxK = 50;
params.missing = -100;

%% INFERENCE
hidden = IBPsampler_run(data,[],params); % no need to initialize Z

%% Complete MATRIX + visualize random image

X_map = IBPsampler_MAP(data.C, hidden.Z, hidden);
% visualization random image
idx = randi(N,1);
figure;
subplot(1,3,1); imagesc(reshape(Xtrue(idx,:),sqrt(784),sqrt(784)) );
subplot(1,3,2); imagesc(reshape(data.X(idx,:),sqrt(784),sqrt(784)) );
subplot(1,3,3); imagesc(reshape(X_map(idx,:),sqrt(784),sqrt(784)) );
%colormap(gray)

disp('SUCCESSFUL');
