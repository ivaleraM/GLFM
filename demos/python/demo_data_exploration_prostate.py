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
from IPython import embed

# import libraries for I/O of data
import cPickle

# ---------------------------------------------
# 1. LOAD DATA TO BE EXPLORED
# ---------------------------------------------
print '\n 1. LOAD DATABASE TO EXPLORE\n'

with open('../../datasets/prostate.pk','rb') as f:
    data = cPickle.load(f)

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tSetting optional parameters for the GLFM model...'

[N, D] = data['X'].shape

params = dict()

# pre-transform a subset of variables
idx_transform = [ D-1 ] # we transform the last dimension
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

params['Niter'] = 100  # number of algorithm iterations (for Gibbs sampler)
params['s2u'] = .005    # Auxiliary variance
params['s2B'] = 1       # Variance of the Gaussian prior of the weigting matrices B
params['alpha'] = 1     # Concentration parameter of the IBP
params['maxK'] = 10     # maximum number of latent features for memory allocation
params['bias'] = 1      # 1 = fix first feature to be active for all patients 

print '\tInitializing Z...'

Kinit = 2   # initial number of latent features
prob = 0.2  # probability of feature activation in matrix Z
hidden = dict()
if params['bias']:
    hidden['Z'] = np.concatenate((np.ones((N,1)),(np.random.rand(N,Kinit-1) < 0.2)*1.0),axis=1)
else:
    hidden['Z'] = (np.random.rand(N,Kinit) < prob) * 1.0

#hidden['Z'] = np.ascontiguousarray( ((np.random.rand(N,Kinit) < prob) * 1.0).astype('float64') )

# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 3. INFERENCE\n'

print '\tInfering latent features...'
tic = time.time()
hidden = GLFM.infer(data,hidden,params=params)
toc = time.time()
time = tic - toc
print '\tElapsed: #.2f seconds.' # (toc-tic)

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
print '\n 4. PROCESSING RESULTS\n'

## Predict MAP estimate for the whole matrix X
patterns = hidden['Z']
X_map = GLFM.computeMAP(data['C'], patterns, hidden, params)

### Compute log lik of the data
loglik = GLFM.compute_log_likelihood(data['X'],data['C'],hidden,params)

## Plot Dimensions
sum(hidden['Z'])
th = 0.03 # threshold to filter out latent features that are not significant
feat_select = np.nonzero(sum(hidden['Z']) >= N*th)[0] # filter features with insufficient number of obs. assigned
hidden['Z']= hidden['Z'][:,feat_select]
hidden['B']= hidden['B'][:,feat_select,:]

sum(hidden['Z'])
[patterns, C, L] = GLFM.get_feature_patterns_sorted(hidden['Z'])

Kest = hidden['B'].shape[1] # number of inferred latent features
D = hidden['B'].shape[0]    # number of dimensions

# choose patterns corresponding to activation of each feature
Zp = np.eye(Kest)
Zp[:,0] = 1 # bias active
leg = ['F0','F1', 'F2', 'F3', 'F4'];
GLFM.plotPatterns(data, hidden, params, Zp, [], [], leg)

print "SUCCESSFUL"

