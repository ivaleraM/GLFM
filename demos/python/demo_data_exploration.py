# ----------------------------------------------------------------
# Demo script for data exploration
# ----------------------------------------------------------------

import numpy as np # import numpy matrix for calculus with matrices
import sys
sys.path.append('../../src/GLFMpython/')
import GLFM        # import General Latent Feature Model Library
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed
import pdb

# import libraries for I/O of data
import scipy.io
import csv

from aux import preprocess
from aux import plot_dim
from aux import plot_dim_1feat

# ---------------------------------------------
# 1. LOAD DATA TO BE EXPLORED
# ---------------------------------------------
print '\n 1. LOAD DATABASE TO EXPLORE\n'

# Fields inside structure
#       xlabel: [502x1 double]
#            X: [502x16 double]
#            C: 'cpncnnccnncpnnpc'
#   cat_labels: {[]  []  []  {10x1 cell}  []  [] {4x1 cell}  ... []}
#       ylabel: {'stage'  'rx'  'dtime' 'status'  'age'  'wt'  'pf'... 'bm'}
#  ylabel_long: {1x16 cell}
input_file = '../../datasets/mat/prostate_v3.mat'
tmp = scipy.io.loadmat(input_file)
data = tmp['data'][0,0] # data is a dictionary with the following keys
(N,D) = data['X'].shape
#X = data['X'].transpose() #  ndarray of dimensions D * N
data['C'] = str(data['C'][0])
pdb.set_trace()

# # dealing with missing data: replace np.nan by -1
# (xx,yy) = np.where(np.isnan(X)) # find positions where X is nan (i.e. missing data)
# for r in xrange(len(xx)):
#     X[xx[r],yy[r]] = -1

idx_toKeep = [1, 2, 4, 13, 15]
data['X'] = data['X'][:,idx_toKeep]
data['C'] = data['C'][idx_toKeep]
data['cat_labels'] = data['cat_labels'][idx_toKeep]
data['ylabel'] = data['ylabel'][idx_toKeep]
data['ylabel_long'] = data['ylabel_long'][idx_toKeep]


# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tInitializing Z...'
Kinit = 1   # initial number of latent features
prob = 0.2  # probability of feature activation in matrix Z
Z = np.ascontiguousarray( ((np.random.rand(Kinit,N) < prob) * 1.0).astype('float64') )
bias = 0
# with bias
#Z = np.concatenate((np.ones((N,1)),(np.random.rand(N,Kest-1) < 0.2)*1.0),axis=1)
#bias = 1

print '\tSetting optional parameters for the GLFM model...'

params['Niter'] = 100  # number of algorithm iterations
params['s2y'] = 0.5    # noise variance for pseudo-obervations
params['s2B'] = 1      # noise variance for feature values
params['s2u'] = 0.1    # auxiliary noise
params['alpha'] = 1    # mass parameter for the Indian Buffet Process

[N, D] = data['X'].shape

# pre-transform a subset of variables
idx_transform = D # we transform the last dimension
params.t = cell(1,D)
params.t_1 = cell(1,D)
params.dt_1 = cell(1,D)
params.ext_dataType = cell(1,D)
for r=idx_transform
    params.t_1{r} = @(x) log(x + 1) # transformation to apply to raw data
    params.t{r} = @(y) exp(y) - 1   # inverse transform to recover raw data
    params.dt_1{r} = @(x) 1./(x+1)  # derivative of inverse transform
    params.ext_dataType{r} = 'p'    # change type of data due to transformation
end

params.missing = -1
params.s2Y = 0         # Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
params.s2u = .005      # Auxiliary variance
if ~isfield(params,'s2B')
    params.s2B = 1     # Variance of the Gaussian prior of the weigting matrices B
end
params.alpha = 1       # Concentration parameter of the IBP
if ~isfield(params,'Niter')
    params.Niter = 1000 # Number of iterations for the gibbs sampler
end
params.maxK = 10
if ~isfield(params,'bias')
    params.bias = 1    # 1 = fix first feature to be active for all patients 
end
params.func = ones(1,D)

if ~isfield(params,'simId')
    params.simId = 1  # simulation identifier (to be able to run the sample sim. multiple times 
end
if ~isfield(params,'save')
    params.save = 0   # save .mat if active
end

## Initialize Hidden Structure

if params.bias
    Zini = [ones(N,1), double(rand(N,1)>0.8)]
else
    Zini = double(rand(N,1)>0.8)
end
hidden['Z'] = Zini # N*K matrix of feature assignments

# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 3. INFERENCE\n'

print '\tInfering latent features...'
tic = time.time()
hidden = GLFM.infer(data,hidden,params)
toc = time.time()
time = tic - toc
print '\tElapsed: #.2f seconds.' # (toc-tic)

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
print '\n 4. PROCESSING RESULTS\n'

Kest = B_out.shape[1] # number of inferred latent features
D = B_out.shape[0]    # number of dimensions


for d in xrange(D):
    ylab = str(data['ylabel_long'][0][d].tolist()[0]) # label for dimension d
    V = np.squeeze(data['cat_labels'][0][d]) # labels for categories (empty if not categorical)
    catlab = tuple( map(lambda x: str(x.tolist()[0]),V) ) # transform list of categories into tuple
    Zp = np.zeros((3,Kest)) # dimensions (numPatterns,Kest)
    Zp[0,0] = 1.0
    Zp[1,1] = 1.0
    Zp[2,2] = 1.0
    plot_dim(X, B_out, Theta_out, C,d,Zp,s2y,s2u,\
            xlabel=ylab, catlabel=catlab) # function to plot patterns corresponding to Zp
    pdb.set_trace()

#k = 1
#d = 3
#ylab = str(data['ylabel_long'][0][d].tolist()[0])
#V = np.squeeze(data['cat_labels'][0][d])
#catlab = tuple( map(lambda x: str(x.tolist()[0]),V) )
#plot_dim_1feat(X, B_out, Theta_out, C,d,k,s2y,s2u, xlabel=ylab, catlabel=catlab)

if params.save
    output_file = [savePath, sprintf('prostateRed_bias#d_simId#d_Niter#d_s2Y#.2f_s2B#.2f_alpha#.2f.mat', ...
        params.bias, params.simId, params.Niter, params.s2Y, params.s2B, params.alpha)]
    save(output_file)
end

## Predict MAP estimate for the whole matrix X
X_map = GLFM_computeMAP(data.C, hidden['Z'], hidden, params)

## Plot Dimensions
if ~params.save

    sum(hidden['Z'])
    th = 0.03 # threshold to filter out latent features that are not significant
    feat_toRemove = find(sum(hidden['Z']) < N*th) # filter features with insufficient number of obs. assigned
    hidden = remove_dims(hidden, feat_toRemove) # remove latent dimensions

    sum(hidden['Z'])
    [patterns, C] = get_feature_patterns(hidden['Z'])

    # choose patterns corresponding to activation of each feature
    Kest = size(hidden.B,2)
    Zp = eye(Kest)
    Zp(:,1) = 1 # bias active
    Zp = Zp(1:min(5,Kest),:)
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5'}

    colors = [] styles = []
    plot_all_dimensions(data, hidden, params, Zp, leg, colors, styles)
end

print "SUCCESSFUL"

