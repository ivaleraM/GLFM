# # DEMO_TOY_IMAGES
## Simple illustration of GLFM pipeline, replicating the example of the
## IBP linear-Gaussian model in (Griffiths and Ghahramani, 2011).

# Let us first import all the necessary libraries.
import numpy as np # import numpy matrix for calculus with matrices
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed
import sys
sys.path.append('../../src/GLFMpython/')
import GLFM

import pdb

# In[2]:

# ---------------------------------------------
# 1. GENERATIVE MODEL
# ---------------------------------------------
print '\n 1. GENERATIVE MODEL\n'

print '\tGenerating feature images (components)...'
# Btrue contains the features images or components in order to generate the
# whole set of images
Btrue = np.array([[0,1.0,0,0,0,0,  1,1,1,0,0,0, 0,1,0,0,0,0,         0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0],         [0,0.0,0,1,1,1,  0,0,0,1,0,1, 0,0,0,1,1,1,         0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0],         [0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,         1,0,0,0,0,0, 1,1,0,0,0,0, 1,1,1,0,0,0],         [0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,         0,0,0,1,1,1, 0,0,0,0,1,0, 0,0,0,0,1,0]])
D = Btrue.shape[1] # number of dimensions
K = Btrue.shape[0] # number of binary images

print '\tPlotting feature images (Btrue)...'
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4] # subfigures handlers
for i in xrange(len(Btrue)):
    pixels = Btrue[i].reshape(int(np.sqrt(Btrue.shape[1])),int(np.sqrt(Btrue.shape[1])))
    # Plot each feature component k
    V[i].imshow(pixels, cmap='gray',interpolation='nearest')
    V[i].set_ylim(0,5) # set Y-axis limits
    V[i].set_xlim(0,5) # set X-axis limits
    V[i].set_title('Image %d' % (i+1)) # set subplot title
#plt.ion()  # turn on interactive mode for plotting (so that the script continues)
plt.show() # display images component
#plt.pause(0.0001)

print '\tSetting model parameters (ground truth) and generate database...'
N = 1000           # number of images to be generated
s2x = 0.5          # noise variance for the observations

print '\tGenerating data with N=%d and noise variance s2x=%.2f' % (N,s2x)
# generate matrix Z
Ztrue = 1.0*(np.random.rand(N,K) < 0.2) #np.random.randint(0,2,size=(N,K)).astype('float64')
# Next line generates the toy database
X = np.sqrt(s2x) * np.random.randn(N,D) + np.inner(Ztrue, Btrue.transpose())



# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------

print '\n 2. INITIALIZATION\n'

print '\tInitializing Z...'
hidden = dict()
Kinit = 2 # initial number of latent features
hidden['Z'] = np.random.randint(0,2,size=(N,Kinit)).astype('float64')

print '\tInitialization of variables needed for the GLFM model...'
data = dict()
data['X'] = X
data['C'] = np.tile('g',(1,X.shape[1]))[0].tostring() # datatype vector

# params is optional
params = dict()
params['alpha'] = 2   # concentration parameter for the IBP
params['Niter'] = 100 # number of algorithm iterations
params['maxK'] = 10


# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\tInfering latent features...\n'
hidden = GLFM.infer(data, hidden, params)

# In[5]:

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
Kest = hidden['B'].shape[1] # number of inferred latent features
D = hidden['B'].shape[0]    # number of dimensions

print '\tPrint inferred latent features...'
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
for k in xrange(hidden['B'].shape[1]):
    if k>len(V):
        break;

    # visualize each inferred dimension
    Zp = np.zeros(Kest)
    Zp[k] = 1.0

    #hidden['B'][:,k]
    #pixels = hidden['B'][:,k].reshape((int(np.sqrt(D)),int(np.sqrt(D))))
    Bpred = GLFM.computeMAP(data['C'],Zp, hidden) # MAP prediction for each dim d
    pixels = Bpred.reshape((int(np.sqrt(D)),int(np.sqrt(D))))
    # Plot
    V[k].imshow(pixels, cmap='gray',interpolation='none')
    V[k].set_ylim(0,5)
    V[k].set_xlim(0,5)
    V[k].set_title('Feature %d' % (k+1))
#plt.ion()  # interactive mode for plotting (script continues)
plt.show() # display figure
#plt.pause(0.0001)

print('\n\n# -------------------')
print "# SUCCESSFUL -- The GLFM is able to recover the original images seamlessly."
print('# -------------------')

