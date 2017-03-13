clear
addpath(genpath('../Ccode/'));
addpath('./aux/');

% Define parameter default values for algorithm
bias = 1; s2Y = 1.0; s2B = 1.0; alpha = 1.0; Niter = 5; maxK = 50; missing_val = -100;

% load MNIST dataset
% This DB can be downloaded here: http://yann.lecun.com/exdb/mnist/
images = loadMNISTImages('../databases/train-images-idx3-ubyte');
N = 100; % subset of images to consider
Xtrue = images(:, randperm(size(images,2),N))' + 1; % N*D matrix (each column is one pixel

perc_missing = 0.2;
mask_missings = rand(size(Xtrue)) < perc_missing;
Xmiss = Xtrue;
Xmiss(mask_missings) = missing_val;

C = repmat('n',1,size(Xmiss,2));

% visualization
idx = randi(N,1);
figure(1); imagesc(reshape(Xmiss(idx,:),sqrt(784),sqrt(784)) );
%colormap(gray)

Xcompl = matrix_completion(Xmiss, C, s2Y, s2B, alpha, Niter, maxK, missing_val);

disp('SUCCESSFUL');
