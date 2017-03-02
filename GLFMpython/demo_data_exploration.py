# ----------------------------------------------------------------
# Demo script for data exploration
# ----------------------------------------------------------------

import numpy as np # import numpy matrix for calculus with matrices
import sys         # import sys to add path of the python wrapper functions
sys.path.append('../Ccode/wrapper_python/')
import GLFM        # import General Latent Feature Model Library
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed
import cPickle     # module to be able to load databases

# import libraries for I/O of data
import scipy.io
import csv

from aux import preprocess
from aux import plot_dim
import pdb

# ---------------------------------------------
# 1. LOAD DATA TO BE EXPLORED
# ---------------------------------------------
print '\n 1. LOAD DATABASE TO EXPLORE\n'

#input_file = '../databases/dataExploration/csv_xls/prostate.csv'
#with open(input_file) as f:
#    reader = csv.reader(f, delimiter=',')
#    for row in reader:
#        print row
#        pdb.set_trace()
#    [Y,genetic_ids, clinical_ids,vocab] = cPickle.load(f)

# Fields inside structure
#       xlabel: [502x1 double]
#            X: [502x16 double]
#            C: 'cpncnnccnncpnnpc'
#   cat_labels: {[]  []  []  {10x1 cell}  []  [] {4x1 cell}  ... []}
#       ylabel: {'stage'  'rx'  'dtime' 'status'  'age'  'wt'  'pf'... 'bm'}
#  ylabel_long: {1x16 cell}
input_file = '../databases/dataExploration/mat/prostate.mat'
tmp = scipy.io.loadmat(input_file)
data = tmp['data'][0,0] # data is a dictionary with the following keys
(N,D) = data['X'].shape
X = data['X'].transpose() #  ndarray of dimensions D * N
C = str(data['C'][0])
# dealing with missing data: replace np.nan by -1
(xx,yy) = np.where(np.isnan(X)) # find positions where X is nan
for r in xrange(len(xx)):
    X[xx[r],yy[r]] = -1

# prepare input data for C++ inference routine # TODO: hide from user
X = preprocess(X,C)

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tInitializing Z...'
Kinit = 1 # initial number of latent features
Z = np.ascontiguousarray( np.random.randint(0,2,size=(Kinit,N)).astype('float64') )

print '\tInitialization of variables needed for the GLFM model...'
# Generate weights for transformation
W = np.ascontiguousarray( 2.0 / np.max(X,1) ) # TODO: account for missings

Niter = 100  # number of algorithm iterations
s2y = 0.5    # noise variance for pseudo-obervations
s2B = 1      # noise variance for feature values
s2u = 0.1    # auxiliary noise
alpha = 1    # mass parameter for the Indian Buffet Process

# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 3. INFERENCE\n'

print '\tInfering latent features...'
tic = time.time()
(Z_out,B_out,Theta_out) = GLFM.infer(X,C,Z,W,Nsim=Niter,s2Y=s2y, s2B=s2B, maxK=D)
toc = time.time()
time = tic - toc
print '\tElapsed: %.2f seconds.' % (toc-tic)

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
print '\n 4. PROCESSING RESULTS\n'

Kest = B_out.shape[1] # number of inferred latent features
D = B_out.shape[0]    # number of dimensions

k = 0
for d in xrange(D):
    # Signature: plot_dim(X,B,Theta,C,d,k,s2Y,s2u,missing=-1,labels=[])
    ylab = str(data['ylabel'][0][d].tolist()[0])
    V = np.squeeze(data['cat_labels'][0][d])
    catlab = tuple( map(lambda x: str(x.tolist()[0]),V) )
    plot_dim(X, B_out, Theta_out, C,d,k,s2y,s2u,\
            xlabel=ylab, catlabel=catlab)
    pdb.set_trace()

print '\tPrint inferred latent features...'
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = \
        plt.subplots(3, 3, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
for k in xrange(B_out.shape[1]):
    if k>len(V):
        break;
    # visualize each inferred dimension
    B_out[:,k]
    pixels = B_out[:,k].reshape((int(np.sqrt(D)),int(np.sqrt(D))))
    pixels
    # Plot
    V[k].imshow(pixels, cmap='gray',interpolation='none')
    V[k].set_ylim(0,5)
    V[k].set_xlim(0,5)
    V[k].set_title('Feature %d' % (k+1))
plt.ion()  # interactive mode for plotting (script continues)
plt.show() # display figure

print "SUCCESSFUL"

