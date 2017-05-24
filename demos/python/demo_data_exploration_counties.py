## --------------------------------------------------
# DEMO: Data exploration on counties database
## --------------------------------------------------

import numpy as np # import numpy matrix for calculus with matrices
import sys
sys.path.append('../../src/GLFMpython/')
import GLFM        # import General Latent Feature Model Library
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed

# import libraries for I/O of data
import cPickle

import pdb

# ---------------------------------------------
# 1. LOAD DATA TO BE EXPLORED
# ---------------------------------------------
print '\n 1. LOAD DATABASE TO EXPLORE\n'

with open('../../datasets/py/counties.pk','rb') as f:
    data = cPickle.load(f)

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tSetting optional parameters for the GLFM model...'

params = dict()
params['Niter'] = 10  # number of algorithm iterations (for Gibbs sampler)
params['s2B'] = 1      # noise variance for feature values
params['s2u'] = 0.005  # auxiliary noise
params['alpha'] = 1    # mass parameter for the Indian Buffet Process

[N, D] = data['X'].shape

# pre-transform a subset of variables
idx_transform = [ 1, 4, 5, 10] # dimensions to be transformed
params['t'] = [None] * D
params['t_1'] = [None] * D
params['dt_1'] = [None] * D
params['ext_dataType'] = [None] * D
for rr in xrange(len(idx_transform)):
    r = idx_transform[rr]
    params['t_1'][r] = lambda x: np.log(x + 1) # transformation to apply to raw data
    params['t'][r] = lambda y: np.exp(y) - 1   # inverse transform to recover raw data
    params['dt_1'][r] = lambda x: 1/(x+1)      # derivative of inverse transform
    params['ext_dataType'][r] = 'p'    # change type of data due to transformation
# dimension 'White' need an inversion too
r = 9
params['t_1'][r] = lambda x: np.log((100-x)+1)  # transformation to apply to raw data
params['t'][r] = lambda y: - np.exp(y) + 101    # inverse transform to recover raw data
params['dt_1'][r] = lambda x: -1/(101 - x)      # derivative of inverse transform
params['ext_dataType'][r] = 'p'                 # change type of data due to transformation

params['maxK'] = 10     # maximum number of latent features for memory allocation
params['bias'] = 1      # 1 = fix first feature to be active for all patients 

print '\tInitializing Z...'

Kinit = 2   # initial number of latent features
prob = 0.2  # probability of feature activation in matrix Z
hidden = dict()
if params['bias']:
    Zini = np.concatenate((np.ones((N,1)),(np.random.rand(N,Kinit-1) < 0.2)*1.0),axis=1)
    #Zini = [ones(N,1), double(rand(N,1)>0.8)]
else:
    #Zini = double(rand(N,1)>0.8)
    hidden['Z'] = (np.random.rand(N,Kinit) < prob) * 1.0
hidden['Z'] = Zini # N*K matrix of feature assignments

#hidden['Z'] = np.ascontiguousarray( ((np.random.rand(N,Kinit) < prob) * 1.0).astype('float64') )

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

Kest = hidden['B'].shape[1] # number of inferred latent features
D = hidden['B'].shape[0]    # number of dimensions

## Predict MAP estimate for the whole matrix X
patterns = hidden['Z']
X_map = GLFM.computeMAP(data['C'], patterns, hidden, params)

## Plot Dimensions
sum(hidden['Z'])
th = 0.03 # threshold to filter out latent features that are not significant
feat_select = np.nonzero(sum(hidden['Z']) >= N*th)[0] # filter features with insufficient number of obs. assigned
hidden['Z']= hidden['Z'][:,feat_select]
hidden['B']= hidden['B'][:,feat_select,:]

sum(hidden['Z'])
[patterns, C, L] = GLFM.get_feature_patterns_sorted(hidden['Z']) # returns sorted patterns

# choose patterns corresponding to activation of each feature
Zp = np.eye(Kest)
Zp[:,1] = 1 # bias active


Zp = patterns[L > 240,:]

colors = np.array([[ 0, 102, 255], [153, 51, 255], \
    [204, 204, 0], [255, 102, 102], \
    [0, 204, 102], [255, 51, 255]])
colors = colors / 255.0
colors[3,:] = [0.9290, 0.6940, 0.1250]
colors[5,:] = [0.4660, 0.6740, 0.1880]

# # change order of colors
# colors = colors([3 5 4 2 1],:)
# colors(4,:) = [0 255 255] ./255

idxD = range(1,D) # idxD = 2:size(data.X,2)
styles = None
#GLFM.plotPatterns(data, hidden, params, Zp)
GLFM.plotPatterns(data, hidden, params, Zp, colors=colors[:Zp.shape[0],:],idxD=idxD)
# plot_all_dimensions(data, hidden, params, Zp, leg, colors, styles)

print "SUCCESSFUL"

