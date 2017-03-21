import numpy as np # package to work with arrays
import random # package used to initialize latent variables
import os
import sys
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])
sys.path.append(root+'/Ccode/wrapper_python/')
import GLFMlib # python wrapper library in order to run C++ inference routine

class GLFM:
    """
    General Latent Feature Model Class
    attributes = model hyperparameters + parameter settings
        alpha : mass parameter for the Indian Buffet Process
        bias  : Number of columns that should be considered as bias (and
                thus, not sampled during inference)
        s2y   : variance of the Gaussian prior for the pseudo-observations
        s2u   : variance of the Gaussian prior on the auxiliary variables
        s2B   : variance of the Gaussian prior on the elements of the
                weighting matrices (latent features) B (TODO: Give intuition)
        Niter : number of internal iterations for the Gibbs sampler within
                the C code before return
        maxK  : max number of latent features to be allocated inside C++
        missing : integer value that should be understood as missing value
        verbose : parameter to control how much info should be printed
    """
    def __init__(self, alpha=1.0, bias=0, s2y=1.0, s2u=0.001, s2B=1.0,\
            Niter=100, maxK=50, missing=-1, verbose=0):
        self.alpha = alpha
        self.bias = bias
        self.s2y = s2y
        self.s2u = s2u
        self.s2B = s2B
        self.Niter = Niter
        self.maxK = maxK
        self.missing = missing
        self.verbose = verbose

    def print_info(self):
        print 'General Latent Feature Model Object:\n'+\
                ('alpha: %.2f\n' % self.alpha)+\
                (' bias: %d\n' % self.bias)+\
                ('  s2y: %.2f\n' % self.s2y)+\
                ('  s2u: %.2f\n' % self.s2u)+\
                ('  s2B: %.2f\n' % self.s2B)+\
                ('Niter: %d\n' % self.Niter)+\
                (' maxK: %d\n' % self.maxK)+\
                ('missing: %d\n' % self.missing)+\
                ('verbose: %d\n' % self.verbose)


    def infer(self, Xin, Cin, Zin):
        # Generate weights for transformation
        Win = np.ascontiguousarray( 2.0 / np.max(Xin,1) )
        return GLFMlib.infer(Xin, Cin, Zin, Win, self.bias, self.s2y, self.s2u,\
            self.s2B, self.alpha, self.Niter, self.maxK, self.missing, self.verbose)

    def complete_matrix(self, Xmiss, C):
        """
        Function to complete missing values of a certain numpy 2dim array

        Input parameters:
             self  : GLFM Object
             Xmiss : numpy array which should be completed, should be normalized (TODO: check)
                     Size [NxD] where N is the number of observations and D is the
                     number of dimensions. Here missing data should be introduced
                     as the numeric value indicated in "missing", and categorical
                     and ordinal data should take values in {1,2, ..., R}.
             C     : char array [1xD] specifying the input data type of each column
                     (dimension) of the observation matrix X. Here 'g' indicates
                     real variable, 'p' positive real variable, 'n' count data,
                     'o' ordinal data and 'c' categorical data.
        Output paramaters:
            Xcompl : same numpy array as Xmiss but whose missing values have been
                     inferred and completed by the algorithm.
        """

        N = Xmiss.shape[0]
        D = Xmiss.shape[1]
        Xmiss[np.isnan(Xmiss)] = missing
        maxK=50 # maximum number of latent features for space allocation

        ## Inference
        #Zini= 1.0*( np.random.rand(N,2) > 0.8 )
        Kinit = 3
        Zini = np.ascontiguousarray( (np.random.rand(Kinit,N) > 0.8).astype('float64') )
        (Zest, B, Theta)= self.infer(Xmiss,C,Zini) # run inference function

        Xcompl=np.copy(Xmiss)
        [idxs_d, idxs_n] = (Xmiss == missing).nonzero()
        f_1= lambda x,w: np.log(np.exp(w*x)-1) # called element-wise 
        f= lambda y,w:  np.log(np.exp(y)+1)/w
        W = 2 / Xmiss.max(1) # D*1 # TODO: Take care of missing values when taking max
        for ii in xrange(len(idxs_n)):
            if Xmiss[idxs_d[ii],idxs_n[ii]] == missing: # will always be the case
                d = idxs_d[ii] # np.ceil(ii/N)
                n = idxs_n[ii] # np.mod(ii,N)
                Br=np.squeeze(B[d,:])
                aux = Zest[:,n].reshape(-1,1) # Zest(:,n)'
                if (C[d] != 'c'):
                    M = np.inner(aux.transpose(),Br)
                if (C[d] == 'g'):
                    # branch checked, working
                    Xcompl[d,n] = M
                elif (C[d] == 'p'):
                    Xcompl[d,n] = f(M,W[d])
                elif (C[d] == 'c'):
                   Br = np.squeeze(B[d,:,:])
                   prob = np.zeros((1,R[d]))
                   Y = np.zeros((1,R[d]))
                   for r in xrange(R[d]):
                       Y[r]= np.inner(aux.transpose(),Br[:,r]) # TODO: check
                   Xcompl[d,n] = np.where(Y == np.max(Y))[0]
                elif (C[d] == 'o' ):
                    #Y = aux * Br # TODO: check dimensions
                    [idx_x, idx_y] = (Theta[d,1:R[d]]>=M).nonzero()
                    Xcompl[d,n] = idx_x[0] # TODO: verify (x,y) and what if more els?
                elif (C[d] == 'n'):
                    Xcompl[d,n] = np.floor(f(M,W[d]))

        return Xcompl

