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
tmp = tmp['data'][0,0]

data = dict() # data is a dictionary with the following keys
data['X'] = tmp[1]
data['C'] = tmp[2]
data['C'] = str(data['C'][0])
data['ylabel'] = tmp[5]
data['cat_labels'] = tmp[3]

idx_toKeep = [0, 1, 3, 12, 14]
data['X'] = data['X'][:,idx_toKeep]
data['C'] = ''.join([data['C'][i] for i in idx_toKeep])
data['ylabel'] = ['Stage','Drug level', 'Prognostic Status',\
        'Size of Primary Tumor (cm^2)', 'Serum Prostatic Acid Phosphatase']
data['cat_labels'] = [['3','4'],['0','1','5'],['alive','vascular',\
        'prostatic ca','others'],None,None]

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tSetting optional parameters for the GLFM model...'

params = dict()
params['Niter'] = 100  # number of algorithm iterations (for Gibbs sampler)
params['s2B'] = 1      # noise variance for feature values
params['s2u'] = 0.1    # auxiliary noise
params['alpha'] = 1    # mass parameter for the Indian Buffet Process

[N, D] = data['X'].shape

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

#params['missing'] = -1  # missing value
params['s2u'] = .005    # Auxiliary variance
params['s2B'] = 1       # Variance of the Gaussian prior of the weigting matrices B
params['alpha'] = 1     # Concentration parameter of the IBP
params['maxK'] = 10     # maximum number of latent features for memory allocation
params['bias'] = 1      # 1 = fix first feature to be active for all patients 
#params['simId'] = 1  # simulation identifier (to be able to run the sample sim. multiple times 

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
hidden = GLFM.infer(data,hidden,params=params)
toc = time.time()
time = tic - toc
print '\tElapsed: #.2f seconds.' # (toc-tic)

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
print '\n 4. PROCESSING RESULTS\n'

Kest = hidden['B'].shape[1] # number of inferred latent features
D = hidden['B'].shape[0]    # number of dimensions

#for d in xrange(D):
#    ylab = str(data['ylabel_long'][0][d].tolist()[0]) # label for dimension d
#    V = np.squeeze(data['cat_labels'][0][d]) # labels for categories (empty if not categorical)
#    catlab = tuple( map(lambda x: str(x.tolist()[0]),V) ) # transform list of categories into tuple
#    Zp = np.zeros((3,Kest)) # dimensions (numPatterns,Kest)
#    Zp[0,0] = 1.0
#    Zp[1,1] = 1.0
#    Zp[2,2] = 1.0
#    plot_dim(X, hidden['B'], hidden['Theta'], C,d,Zp,s2y,s2u,\
#            xlabel=ylab, catlabel=catlab) # function to plot patterns corresponding to Zp
#    pdb.set_trace()
#
#k = 1
#d = 3
#ylab = str(data['ylabel_long'][0][d].tolist()[0])
#V = np.squeeze(data['cat_labels'][0][d])
#catlab = tuple( map(lambda x: str(x.tolist()[0]),V) )
#plot_dim_1feat(X, hidden['B'], hidden['Theta'], C,d,k,s2y,s2u, xlabel=ylab, catlabel=catlab)

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
[patterns, C, L] = GLFM.get_feature_patterns(hidden['Z'])

# choose patterns corresponding to activation of each feature
Zp = np.eye(Kest)
Zp[:,1] = 1 # bias active
Zp = Zp[1:min(5,Kest),:]

colors = None
styles = None
pdb.set_trace()
# plot_all_dimensions(data, hidden, params, Zp, leg, colors, styles)

print "SUCCESSFUL"

