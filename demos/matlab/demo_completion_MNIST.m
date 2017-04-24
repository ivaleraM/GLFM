%% demo TOY COMPLETE MNIST

clear
addpath(genpath('../../src/')); % location of GLFM source code
rng( round(sum(1e5*clock)) );

%% LOAD MNIST dataset
% This DB can be downloaded here: http://yann.lecun.com/exdb/mnist/
images = loadMNISTImages('../../datasets/raw/MNIST/train-images-idx3-ubyte');
N = 50;                                             % subset of images to consider
% Xtrue is the ground truth data (each column is one pixel)
Xtrue = images(:, randperm(size(images,2),N))' + 1; % N*D matrix

% Build Xmiss (data with missing values)
perc_missing = 0.2;
missing_val = -100;
mask_missings = rand(size(Xtrue)) < perc_missing;
Xmiss = Xtrue;
Xmiss(mask_missings) = missing_val;

data.X = Xmiss;
data.C = repmat('n',1,size(Xmiss,2)); % data type, here pixel values onsidered as count data

%% Define parameter default values for algorithm
params.bias = 1;        % bias feature (active for all obs. to capture mean effects)
%params.s2Y = 0;         % variance noise for pseudo-observations
params.s2B = 1;         % variance noise for feature matrix
params.alpha = 10;      % concentration parameter for the Indian Buffet Process
params.Niter = 10;      % number of iterations to train the GLFM model
params.maxK = 50;       % maximum number of latent features (for memory allocation)
params.missing = -100;  % values for missings

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

disp('SUCCESSFUL');
