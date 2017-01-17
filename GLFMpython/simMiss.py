# example python script
import numpy as np
import time
import scipy.io
from matrix_completion import matrix_completion

missing = -1
p = 0.1 # probability of missing
f_1 = lambda x,w: np.log(np.exp(w*x)-1)
df_1 = lambda x,w: w/(1-np.exp(-w*x))
logN = lambda x: np.max(np.log(x),-30) # To avoid -Inf

## Selecting missing data
np.random.seed( int(time.time()) )
print 'loading Wine.mat database...'
tmp = scipy.io.loadmat('../databases/Wine.mat')
C = tmp['C']
D = tmp['D']
N = tmp['N']
R = tmp['R']
W = tmp['W']
X = tmp['X']
maxR = tmp['maxR']
del tmp

Nmiss = np.round(N*D*p)
miss = np.random.permutation(N*D)
miss = miss[:Nmiss]
Xmiss = X # Observation matrix
X_tmp = Xmiss.flatten()
X_tmp[miss] = missing # Missing data are coded as missing
Xmiss = X_tmp.reshape(Xmiss.shape[0],-1)
del X_tmp
Xmiss[np.isnan(X)] = missing
N = X.shape[0]
D = X.shape[1]
s2Y = 1   # Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
s2B = 1   # Variance of the Gaussian prior of the weigting matrices B
alpha = 1 # Concentration parameter of the IBP
Nsim = 5000 # Number of iterations for the gibbs sampler
bias = 0
maxK = D
## Inference
Zini = 1.0*( np.random.rand(N,2) > 0.8 )
# Call inner C function
[Zest B Theta]= IBPsampler(Xmiss,C,Zini,bias,s2Y,s2B,alpha,Nsim,maxK,missing)


## Compute test log-likelihood
# TODO: Initialize structure TLK
XT = X
idxs = (Xmiss == missing).nonzero()
ii = 0
for i in xrange(idxs[0].shape[0]): # for each missing element
    ii=ii+1
    if (XT[n,d]~=-1 || ~isnan(XT[n,d])) # TODO: Why do we need this condition?
        n = idxs[0][i]
        d = idxs[1][i]

        j = j+1
        Br = np.squeeze(B[d,:,0])
        aux = Zest[:,n].reshape(-1,1) # aux 
        if (C[d] == 'g'):
            TLK[n,d] = logN(normpdf(XT[n,d],aux * Br,1))

        elif (C[d] == 'p'):
            TLK[n,d] = logN(normpdf(f_1(XT[n,d],W[d]),aux * Br,2)) + \
                    logN(np.abs(df_1(XT[n,d],W[d])))

        elif (C[d] == 'c'):
            Br = np.squeeze(B[d,:,:])
            prob = np.zeros((1,R[d]))
            for r in xrange(R[d]-1):
                for r2 in xrange(R[d]):
                    if r2!=r:
                        prob[r] =prob[r]+ logN(normcdf(aux *(Br(:,r)-Br(:,r2)),0,1))
            if (1-np.sum(np.exp(prob[0:R[d]-1])) < 0):
                 prob[-1] = -30
                 prob = logN(np.exp(prob) / np.sum(np.exp(prob))) % TODO: Matrix division!!
            else:
                 prob[-1] = logN( 1-np.sum(np.exp(prob[:(R[d]-1)])) )
            TLK[n,d] =prob[XT[n,d]]

        elif (C[d] == 'o'):
            Br = np.squeeze(B[d,:,0])
            if XT[n,d] == 1:
                TLK[n,d] = logN(normcdf(Theta[d,XT[n,d]]- aux*Br,0,1))
            elif XT[n,d] == R[d]:
                # TODO: Check index of Theta? Is -1 necessary? Below again.
                TLK[n,d] = logN(1- normcdf(Theta[d,XT[n,d]-1] - aux*Br,0,1))
            else:
                TLK[n,d] = logN(normcdf(Theta[d,XT[n,d]] - aux*Br,0,1) \
                        - normcdf(Thetar[d,XT[n,d]-1] - aux*Br,0,1))
        elif (C[d] == 'n'):
            Br = np.squeeze(B[d,:,0])
            if (XT[n,d] == 0): # TODO: Check meaning?
                TLK[n,d] = logN(normcdf(f_1(XT[n,d]+1,W[d]) - aux*Br,0,1))
            else:
                TLK[n,d] = logN(normcdf(f_1(XT[n,d]+1,W[d]) - aux*Br,0,1) \
                        - normcdf(f_1(XT[n,d],W[d]) - aux*Br,0,1))
