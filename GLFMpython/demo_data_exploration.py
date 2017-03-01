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
import scipy.io

# ---------------------------------------------
# 1. LOAD DATA TO BE EXPLORED
# ---------------------------------------------
print '\n 1. LOAD DATABASE TO EXPLORE\n'

input_file = ''
with open(input_file) as f:
    [Y,genetic_ids, clinical_ids,vocab] = cPickle.load(f)

data = scipy.io.loadmat(input_file)

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 2. INITIALIZATION\n'

print '\tInitializing Z...'
Kinit = 1 # initial number of latent features
Z = np.ascontiguousarray( np.random.randint(0,2,size=(Kinit,N)).astype('float64') )

print '\tInitialization of variables needed for the GLFM model...'
C = np.tile('g',(1,X.shape[0]))[0].tostring() # vector to indicate datatype of each dimension
# Generate weights for transformation
W = np.ascontiguousarray( 2.0 / np.max(X,1) )

Niter = 100  # number of algorithm iterations
s2y = s2x/5  # noise variance for pseudo-obervations
s2B = 1      # noise variance for feature values
s2u = 0.1    # auxiliary noise
alpha = 1    # mass parameter for the Indian Buffet Process

# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\n 3. INFERENCE\n'

print '\tInfering latent features...'
tic = time.time()
(Z_out,B_out,Theta_out) = GLFM.infer(X,C,Z,W,Nsim=Niter,s2Y=s2y, s2B=s2B, maxK=10)
toc = time.time()
time = tic - toc
print '\tElapsed: %.2f seconds.' % (toc-tic)

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
print '\n 4. PROCESSING RESULTS\n'

Kest = B_out.shape[1] # number of inferred latent features
D = B_out.shape[0]    # number of dimensions

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

