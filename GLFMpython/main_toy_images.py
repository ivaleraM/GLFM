import numpy as np
import sys
sys.path.append('../Ccode/wrapper_python/')
import GLFM
import matplotlib.pyplot as plt
import time

import pdb

print 'Generating toy_images database...'
Btrue = np.array([[0,1.0,0,0,0,0,  1,1,1,0,0,0, 0,1,0,0,0,0, \
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0], \
        [0,0.0,0,1,1,1,  0,0,0,1,0,1, 0,0,0,1,1,1, \
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0], \
        [0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, \
        1,0,0,0,0,0, 1,1,0,0,0,0, 1,1,1,0,0,0], \
        [0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, \
        0,0,0,1,1,1, 0,0,0,0,1,0, 0,0,0,0,1,0]])
Btrue = Btrue

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4]
for i in xrange(len(Btrue)):
    pixels = Btrue[i].reshape((np.sqrt(Btrue.shape[1]),np.sqrt(Btrue.shape[1])))
    # Plot
    V[i].imshow(pixels, cmap='gray',interpolation='none')
    V[i].set_ylim(0,5)
    V[i].set_xlim(0,5)
    V[i].set_title('Image %d' % (i+1))
plt.ion() # interactive mode for plotting (script continues)
plt.show()

print 'Initialization and parameter settings...'
N = 500 # number of images to be generated
D = Btrue.shape[1]
K = Btrue.shape[0] # number of Binary images
Ztrue = np.ascontiguousarray( np.random.randint(0,2,size=(K,N)).astype('float64') )
s2x = 0.5

X = np.sqrt(s2x) * np.random.randn(D,N) + np.inner(Btrue.transpose(),Ztrue.transpose())
X = np.ascontiguousarray(X)
#X = X - 0.5 # center data
C = np.tile('g',(1,X.shape[0]))[0].tostring()
Kinit = 1
Z = np.ascontiguousarray( np.random.randint(0,2,size=(Kinit,N)).astype('float64') )

print 'Infering latent features...'
tic = time.time()
(Z_out,B_out,Theta_out) = GLFM.infer(X,C,Z,Nsim=1000,s2Y=s2x/5, s2B=1, maxK=10)
toc = time.time()
time = tic - toc
print 'Elapsed: %.2f seconds' % (toc-tic)

Kest = B_out.shape[1]
D = B_out.shape[0]

f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = \
        plt.subplots(3, 3, sharex='col', sharey='row')
V = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
for k in xrange(B_out.shape[1]):
    if k>len(V):
        break;
    # visualize each inferred dimension
    B_out[:,k]
    pixels = B_out[:,k].reshape((np.sqrt(D),np.sqrt(D)))
    pixels
    # Plot
    V[k].imshow(pixels, cmap='gray',interpolation='none')
    V[k].set_ylim(0,5)
    V[k].set_xlim(0,5)
    V[k].set_title('Feature %d' % (k+1))
plt.ion() # interactive mode for plotting (script continues)
plt.show()

print "SUCCESSFUL"

