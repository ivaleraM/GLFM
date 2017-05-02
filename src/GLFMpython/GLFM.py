import numpy as np
import random

import mapping_functions as mf
import matplotlib.pyplot as plt

import time as timeI
import os
import sys
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])
sys.path.append(os.path.join(root, 'Ccode/wrapper_python/'))

import GLFMlib # python wrapper library in order to run C++ inference routine
import mapping_functions as mf

def infer(Xin,Cin,Zin,bias=0, s2u=0.001, s2B=1.0,
        alpha=1.0, Nsim=100, maxK=50, missing=-1, verbose=0):
    """
    Python wrapper to launch inference routine for the GLFM model.
    Inputs:
        Xin: input data D*N where:
                N = number of observations
                D = number of dimensions
        Cin: array indicating types of data ('g': real,'p': positive real-valued,
            'c': categorical 'o': ordinal 'n': count data)
        Zin: initial feature activation matrix: K*N
                K = number of latent dimensions
        bias: indicator of whether to include or not a bias
        s2u: internal auxiliary noise
        s2B: noise variance for prior over elements of matrix B
        alpha: concentration parameter of the Indian Buffet Process
        Nsim: number of simulations
        maxK: maximum number of features for memory allocation
        missing: value for missings (should be an integer, not nan)
        verbose: indicator to print more information
    Output:
        Z_out: feature activation matrix sampled from posterior
        B_out: observation matrix sampled from posterior
        Theta_out: auxiliary variables for ordinal data (needed to compute MAP,
                    or posterior PDFs)
        mu_out: mean parameter for internal transformation
        w_out: scale parameter for internal transformation
        s2Y_out: inferred noise variance for pseudo-observations Y
    """
    # prepare input data for C++ inference routine
    Fin = 1.0 * np.ones(Xin.shape[0]) # choose internal transform function (for positive)
    Xin = np.ascontiguousarray( Xin ) # specify way to store matrices to be
    Zin = np.ascontiguousarray( Zin ) # compatible with C code
    tic = timeI.time()
    # RUN C++ routine
    (Z_out,B_out,Theta_out,mu_out,w_out,s2Y_out) = \
            GLFMlib.infer(Xin, Cin, Zin, Fin, bias, s2u, s2B, alpha, Nsim,\
        maxK, missing, verbose)
    toc = timeI.time()
    time = tic - toc
    print '\tElapsed: #.2f seconds.' # (toc-tic)
    return (Z_out,B_out,Theta_out,mu_out,w_out,s2Y_out)

def complete(Xmiss, C, bias=0, s2Y=1, s2u=1, s2B=1, alpha=1, Niter=50, missing=-1):
    """
    Function to complete missing values of a certain numpy 2dim array

    Input parameters:
         Xmiss : numpy array which should be completed.
                 Size [NxD] where N is the number of observations and D is the
                 number of dimensions. Here missing data should be introduced
                 as the numeric value indicated in "missing".
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
    Xmiss[np.isnan(Xmiss)] = missing
    maxK=50 # maximum number of latent features for space allocation

    ## Inference
    #Zini= 1.0*( np.random.rand(N,2) > 0.8 )
    Kinit = 3
    Zini = (np.random.rand(Kinit,N) > 0.8).astype('float64')
    # Call inner C function
    (Zest, B, Theta)= infer(Xmiss,C,Zini,bias,s2u,s2B,alpha,Niter,maxK,missing)

    Xcompl=np.copy(Xmiss)
    [idxs_d, idxs_n] = (Xmiss == missing).nonzero()

    for ii in xrange(len(idxs_n)): # for each missing
        if Xmiss[idxs_d[ii],idxs_n[ii]] == missing: # will always be the case
            d = idxs_d[ii]
            n = idxs_n[ii]
            Br=np.squeeze(B[d,:])
            aux = Zest[:,n].reshape(-1,1) # Zest(:,n)'
            M = np.inner(aux.transpose(),Br)
            if (C[d] != 'c'):
                a = 1
                # TODO
            if (C[d] == 'g'):
                # branch checked, working
                Xcompl[d,n] = M
            elif (C[d] == 'p'):
                Xcompl[d,n] = mf.fpos(M)
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
                Xcompl[d,n] = np.floor(mf.fpos(M))
    return Xcompl

def computeMAP(C, Zp, hidden, params, idxsD=[]):
    """
    Function to generate the MAP solution corresponding to patterns in Zp
    Inputs:
      C: 1*D string with data types, D = number of dimensions
      Zp: P * K matrix of feature activation for which to compute the MAP estimate
          (P is the number of obs.)
      hidden: structure with latent variables learned by the model
          - B: latent feature matrix (D * K * maxR)  where
                  D: number of dimensions
                  K: number of latent variables
               maxR: maximum number of categories across all dimensions
          - mu: 1*D shift parameter
          - w:  1*D scale parameter
          - theta: D*maxR matrix of auxiliary vars (for ordinal variables)
    ----------------(optional) ------------------
          - idxsD: dimensions to infer

    Outputs:
      X_map: P*Di matrix with MAP estimate where Di = length(idxsD)
    """

    if (len(idxsD) == 0):
        idxsD = range(hidden['B'].shape[0])

    P = Zp.shape[0]
    K = hidden['B'].shape[1]
    assert(Zp.shape[1] ~= K), "Incongruent sizes between Zp and hidden.B: number of latent variables should not be different"

    X_map = np.zeros((P,len(idxsD))) # output
    for dd in xrange(len(idxsD)): # for each dimension
        d = idxsD[dd]
        # if params.has_key('t'): # if external transformations have been defined
        #     if len(params['t'][d]) > 0: # there is an external transform for data type d
        #         C(d) = params.ext_dataType{d} # set new type of data

        if C[d] == 'g':
            X_map[:,d] = mf.f_g( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) )
        elif C[d] == 'p':
            X_map(:,d) = f_p( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) )
        elif C[d] == 'n':
            X_map(:,d) = f_n( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) )
        elif C[d] == 'c':
            X_map(:,d) = f_c( Zp * squeeze(hidden.B(d,:,1:hidden.R(d))) )
        elif C[d] == 'o':
            X_map(:,d) = f_o( Zp * squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)) )
        else:
            raise ValueError('Unknown data type')
        end
        if (sum(np.isnan(X_map[:,d])) > 0):
            raise ValueError('Some values are nan!')
        # if params.has_key('t'):
        #     if len(params['t'][d]) > 0: # there is an external transform for data type d
        #         X_map[:,d] = params['t'][d]( X_map[:,d] )
    return X_map

def computePDF(data, Zp, B, s2y, missing_val=-1,  params, d):
    """
    Function to compute probability density function for dimension d
    """
    X[np.isnan(X[:,d]),d] = missing_val

    # compute x-domain [mm MM] to compute pdf
    mm = np.min(X[X[:,d] ~= missing_val, d]) # min value
    MM = np.max(X[X[:,d] ~= params.missing, d]) # max value

    if (len(params['t'][d]) == 0): # if there is an external transformation
        C[d] = params['ext_dataType'][d]
        mm = params['t_1'][d][mm]
        MM = params['t_1'][d][MM]

    P = Zp.shape[0]
    K = B.shape[1]
    assert (Zp.shape[1] == K), "Incongruent sizes between Zp and hidden.B"
    if (C[d] == 'g') || (C[d] == 'p'):
        if not(params.has_key('numS')):
            numS = 100
        xd = np.linspace(mm, MM, num=numS)
    elseif (C[d] == 'n'):
        xd = range(mm,MM+1)
        numS = len(xd)
    else:
        xd = np.unique(X[X[:,d] ~= missing_val, d])
        numS = length(xd) # number of labels for categories or ordinal data
    pdf = np.zeros((P,numS))
    for p in xrange(P):
        if C[d] == 'g':
            pdf[p,:] = mf.pdf_g(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params)
        elif C[d] == 'p':
            pdf[p,:] = mf.pdf_p(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params)
        elif C[d] == 'n':
            pdf[p,:] = mf.pdf_n(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params)
        elif C[d] == 'c':
            pdf[p,:] = mf.pdf_c(Zp(p,:), squeeze(hidden.B(d,:,1:hidden.R(d))), hidden.s2Y(d))
        elif C[d] == 'o':
            pdf[p,:] = mf.pdf_o(Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)), hidden.s2Y(d))
        else:
                error('Unknown data type')
        assert (np.sum(np.isnan(pdf)) == 0), "Some values are nan!"

    if params.has_key('t'):
        if len(params['t'] == 0: # we have used a special transform beforehand
            xd = params['t'][d](xd) # if there was an external transformation, transform pdf
            pdf = pdf .* np.abs( params['dt_1'][d](xd) )



def plot_dim_1feat(X,B,Theta,C,d,k,s2Y,s2u,missing=-1,catlabel=[],xlabel=[]):
    """
    Function to plot an individual dimension of feature matrix B
    Inputs:
        X: observation matrix of dimensions (D,N)
        B: (D,Kest,maxR) ndarray
        C: datatype vector - str of length D
        d: dimension to plot
        k: feature to consider
    Output:
        void
    """
    plt.figure()
    plt.xlabel(xlabel)

    (D,Kest,maxR) = B.shape
    Xd = X[d,:]
    Cd = C[d]
    if k<0 or k>Kest:
        print('Error: k index should be bigger than o and smaller than Kest')
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    if Cd == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
        plt.plot(xx,pdf)

    elif Cd == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w))
        plt.plot(xx,pdf)

    elif Cd == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
        plt.stem(xx,pdf)

    elif Cd == 'c':
        R = len( np.unique(Xd[mask]) )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = np.squeeze(B[d,:,:]) # TODO: Review that size = [K*R]
        pdf = mf.pdf_cat(Zn,Bdv,s2u,R)
        bar_width = 0.35
        index = np.arange(len(pdf))
        plt.bar(index,pdf,bar_width)
        plt.xticks(index + bar_width / 2, catlabel, rotation='vertical')

    elif Cd == 'o':
        a = 1
        # TODO
    else:
        print 'Unknown datatype'
    plt.ion()
    plt.show()
    plt.pause(0.0001)
    return

def plot_dim(X,B,Theta,C,d,Zp,s2Y,s2u,missing=-1,catlabel=[],xlabel=[]):
    """
    Function to plot an individual dimension of feature matrix B
    Inputs:
        X: observation matrix of dimensions (D,N)
        B: (D,Kest,maxR) ndarray
        C: datatype vector - str of length D
        d: dimension to plot
    Output:
        void
    """
    if (Zp.shape[1] != B.shape[1]):
        print 'Error: Sizes of Zp and B are inconsistent'

    colors = ['r','b','g','m','g']
    plt.figure()       # create new figure
    plt.xlabel(xlabel) # add x legend
    #print xlabel

    (D,Kest,maxR) = B.shape
    Xd = X[d,:]
    Cd = C[d]
    # only consider values in dimension d which are not missing
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    (numPatterns,Kest) = Zp.shape
    if Cd == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
            plt.plot(xx,pdf,label=str(Zn))

    elif Cd == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w))
            plt.plot(xx,pdf,colors[p],label=str(Zn))

    elif Cd == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
            plt.stem(xx,pdf,colors[p], label=str(Zn))

    elif Cd == 'c':
        R = len( np.unique(Xd[mask]) )
        Bdv = np.squeeze(B[d,:,:]) # TODO: Review that size = [K*R]
        bar_width = 0.6/numPatterns
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_cat(Zn,Bdv,s2u,R)
            index = np.arange(len(pdf))
            plt.bar(index+p*bar_width,pdf,width=bar_width,color=colors[p]) #, label=str(Zn))
        plt.xticks(index + bar_width / 2, catlabel, rotation='vertical')
#ax.bar(x-0.2, y,width=0.2,color='b',align='center')
#ax.bar(x, z,width=0.2,color='g',align='center')

    elif Cd == 'o':
        print 'This category is currently under development'
        a = 1
        # TODO
    else:
        print 'Unknown datatype'
    plt.legend()
    #plt.ion()
    plt.show()
    plt.pause(0.0001)

    return
