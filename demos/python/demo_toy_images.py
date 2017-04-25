
# coding: utf-8

# # Introduction to the General Latent Feature Model (GLFM)
# # DEMO_TOY_IMAGES

# Let us first import all the necessary libraries.

# In[1]:

import numpy as np # import numpy matrix for calculus with matrices
import matplotlib.pyplot as plt # import plotting library
import time        # import time to be able to measure iteration speed
import sys
sys.path.append('../../src/GLFMpython/')
import GLFM


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
plt.ion()  # turn on interactive mode for plotting (so that the script continues)
plt.show() # display images component
plt.pause(0.0001)

print '\tSetting model parameters (ground truth) and generate database...'
N = 1000           # number of images to be generated
s2x = 1            # noise variance for the observations

print '\tGenerating data with N=%d and noise variance s2x=%.2f' % (N,s2x)
# generate matrix Z
Ztrue = np.random.randint(0,2,size=(K,N)).astype('float64')
# Next line generates the toy database
X = np.sqrt(s2x) * np.random.randn(D,N) + np.inner(Btrue.transpose(),Ztrue.transpose())


# In[3]:

# ---------------------------------------------
# 2. INITIALIZATION FOR GLFM ALGORITHM
# ---------------------------------------------

print '\n 2. INITIALIZATION\n'

print '\tInitializing Z...'
Kinit = 1 # initial number of latent features
Z = np.random.randint(0,2,size=(Kinit,N)).astype('float64')

print '\tInitialization of variables needed for the GLFM model...'
C = np.tile('g',(1,X.shape[0]))[0].tostring() # vector to indicate datatype of each dimension

Niter = 100  # number of algorithm iterations
s2y = 1      # noise variance for pseudo-obervations
s2B = 1      # noise variance for feature values



# In[4]:

# ---------------------------------------------
# 3. RUN INFERENCE FOR GLFM ALGORITHM
# ---------------------------------------------
print '\tInfering latent features...'
(Z_out,B_out,Theta_out) = GLFM.infer(X,C,Z,Nsim=Niter, s2B=s2B, s2u=0.005, maxK=10, transform='off')


# In[5]:

# ---------------------------------------------
# 4. PROCESS RESULTS
# ---------------------------------------------
Kest = B_out.shape[1] # number of inferred latent features
D = B_out.shape[0]    # number of dimensions

print '\tPrint inferred latent features...'
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) =         plt.subplots(3, 3, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
for k in xrange(B_out.shape[1]):
    if k>len(V):
        break;
        
    # visualize each inferred dimension
    Zp = np.zeros(Kest)
    Zp[k] = 1.0
    Bpred = GLFM.predict(X,C,Zp,B_out,Theta_out) # MAP prediction for each dimension d, X necessary to know suffStats
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
plt.pause(0.0001)

print('\n\n# -------------------')
print "# SUCCESSFUL"
print('# -------------------')

