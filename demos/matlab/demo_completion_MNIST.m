%% Demo TOY COMPLETE MNIST
% Illustration of missing data estimation on the MNIST image dataset.
close 
clear
addpath(genpath('../../src/')); % location of GLFM source code
rng( round(sum(1e5*clock)) );

%% LOAD MNIST dataset
% This DB can be downloaded here: http://yann.lecun.com/exdb/mnist/
images = loadMNISTImages('../../datasets/raw/MNIST/train-images-idx3-ubyte');
N = 50;                                             % subset of images to consider
% Xtrue is the ground truth data (each column is one pixel)
Xtrue = images(:, randperm(size(images,2),N))' + 1; % N*D matrix

% Build Xmiss (observation matrix with missing values)
perc_missing = 0.3; % total percentage of missings
missing_val = -100; % value that encodes missing
mask_missings = rand(size(Xtrue)) < perc_missing;
Xmiss = Xtrue;
Xmiss(mask_missings) = missing_val;

data.X = Xmiss;
data.C = repmat('n',1,size(Xmiss,2)); % data type, here pixel values onsidered as count data
% idxOrd = find(max(data.X,[],1)>1) ; %
% data.C(idxOrd)='o';

%% Define parameter default values for algorithm
params.bias = 0;        % bias feature (active for all obs. to capture mean effects)
params.s2B = 0.5;         % variance noise for feature matrix
params.alpha = 10;      % concentration parameter for the Indian Buffet Process
params.Niter = 50;      % number of iterations to train the GLFM model
params.maxK = 100;      % maximum number of latent features (for memory allocation)
params.missing = -100;  % values for missings


%% INFERENCE
[Xcompl, hidden] = GLFM_complete(data,[],params); % no need to initialize Z
% Alternatively, we can compute the whole MAP solution X_map
% hidden = GLFM_infer(data,[],params); % no need to initialize Z
% X_map = GLFM_computeMAP(data.C, hidden.Z, hidden, params);

%% Visualize random image
idx = randi(N,1);
figure;
subplot(1,3,1); imagesc(reshape(Xtrue(idx,:),sqrt(784),sqrt(784)) );
subplot(1,3,2); imagesc(reshape(data.X(idx,:),sqrt(784),sqrt(784)) );
subplot(1,3,3); imagesc(reshape(Xcompl(idx,:),sqrt(784),sqrt(784)) );

disp('SUCCESSFUL');
