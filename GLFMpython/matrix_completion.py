import numpy as np
import random

import sys
sys.path.append('../Ccode/wrapper_python/')
import GLFM

import pdb

def complete_matrix(Xmiss, C, bias=0, s2Y=1, s2B=1, alpha=1, Niter=50, missing=-1):
    """
    Function to complete missing values of a certain numpy 2dim array

    Input parameters:
         Xmiss : numpy array which should be completed, should be normalized (TODO: check)
                 Size [NxD] where N is the number of observations and D is the
                 number of dimensions. Here missing data should be introduced
                 as the numeric value indicated in "missing", and categorical
                 and ordinal data should take values in {1,2, ..., R}.
         C     : char array [1xD] specifying the input data type of each column
                 (dimension) of the observation matrix X. Here 'g' indicates
                 real variable, 'p' positive real variable, 'n' count data,
                 'o' ordinal data and 'c' categorical data.
         bias  : Number of columns that should be considered as bias (and
                 thus, not sampled during inference)
         s2Y   : variance of the Gaussian prior on the auxiliary variables
                 (pseudo-observations) Y (TODO: detail how to change it)
         s2B   : variance of the Gaussian prior on the elements of the
                 weighting matrices (latent features) B (TODO: Give intuition)
         alpha : mass parameter for the Indian Buffet Process
         Niter : number of internal iterations for the Gibbs sampler within
                 the C code before return
         missing : integer value that should be understood as missing value
    Output paramaters:
        Xcompl : same numpy array as Xmiss but whose missing values have been
                 inferred and completed by the algorithm.
    """

    N = Xmiss.shape[0]
    D = Xmiss.shape[1]
    pdb.set_trace()
    Xmiss[np.isnan(Xmiss)] = missing
    maxK=50 # maximum number of latent features for space allocation

    ## Inference
    #Zini= 1.0*( np.random.rand(N,2) > 0.8 )
    Kinit = 3
    Zini = np.ascontiguousarray( (np.random.rand(Kinit,N) > 0.8).astype('float64') )
    # Call inner C function
    (Zest, B, Theta)= GLFM.infer(Xmiss,C,Zini,bias,s2Y,s2B,alpha,Niter,maxK,missing)

    # Compute test log-likelihood
    # TODO: Compute test LLH?
    Xcompl=Xmiss
    idxs = (Xmiss == missing).nonzero()
    #miss=find(Xmiss==missing)';
    f_1= lambda x,w: np.log(np.exp(w*x)-1) # TODO: Verify if called element-wise or not
    f= lambda y,w:  np.log(np.exp(y)+1)/w

    for ii in xrange(len(idxs)):
        if Xmiss[idxs[ii][0],idxs[ii][1]] == missing: # will always be the case
            d = idxs[ii][1] # np.ceil(ii/N)
            n = idxs[ii][0] # np.mod(ii,N)
            #if (n==0)
            #    n=N;
            #end
            Br=np.squeeze(B[d,:,1]) # TODO: Check dimensions of matrix B
            aux = Zest[:,n].reshape(-1,1) # Zest(:,n)'
            if (C[d] == 'g'):
                Xcompl[n,d] = aux * Br
            elif (C[d] == 'p'):
                Xcompl[n,d] = f(aux * Br,W[d])
            elif (C[d] == 'c'):
               Br = np.squeeze(B[d,:,:])
               prob = np.zeros((1,R[d]))
               Y = np.zeros((1,R[d]))
               for r in xrange(R[d]):
                   Y[r]= aux * Br[:,r]
               Xcompl[n,d] = np.where(Y == np.max(Y))[0]
            elif (C[d] == 'o' ):
                Br = np.squeeze(B[d,:,1])
                Y = aux * Br # TODO: check dimensions
                [idx_x, idx_y] = (Theta[d,1:R[d]]>=Y).nonzero()
                Xcompl[n,d] = idx_x[0] # TODO: verify (x,y) and what if more els?
            elif (C[d] == 'n'):
                Br = np.squeeze(B[d,:,1])
                Xcompl[n,d] = np.floor(f(aux * Br,W[d]))

    return Xcompl
