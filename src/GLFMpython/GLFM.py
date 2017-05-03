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

import pdb

def infer(data,hidden=dict(), params=dict()):
    """
    Python wrapper to launch inference routine for the GLFM model.
    Inputs:
        data: dictionary containing all input data structures
            X: input data N*D where:
                N = number of observations
                D = number of dimensions
            C: array indicating types of data ('g': real,'p': positive real-valued,
                'c': categorical; 'o': ordinal; 'n': count data)
        hidden (optional): dictionary containing latent variables
            Z: initial feature activation matrix: N*K
                K = number of latent dimensions
        params (optional): dictionary containing all simulation parameters
            bias: indicator of whether to include or not a bias
            s2u: internal auxiliary noise
            s2B: noise variance for prior over elements of matrix B
            alpha: concentration parameter of the Indian Buffet Process
            Nsim: number of simulations
            maxK: maximum number of features for memory allocation
            missing: value for missings (should be an integer, not nan)
            verbose: indicator to print more information

    Output:
        hidden:
            Z: feature activation matrix sampled from posterior
            B: observation matrix sampled from posterior
            Theta: auxiliary variables for ordinal data (needed to compute MAP,
                    or posterior PDFs)
            mu: mean parameter for internal transformation
            w: scale parameter for internal transformation
            s2Y: inferred noise variance for pseudo-observations Y
    """
    # complete dictionary params
    params = init_default_params(data, params) # complete unspecified fields

    # check input syntax
    assert(type(data) is dict), "input 'data' should be a dictionary."
    assert(type(hidden) is dict), "input 'hidden' should be a dictionary."
    assert(type(params) is dict), "input 'params' should be a dictionary."
    assert(data.has_key('X')), "dictionary data does not have any matrix X defined."
    assert(data.has_key('C')), "dictionary data does not have any datatype array C defined."
    assert(params['bias'] <= 1), "bias parameter misspecified: should be either 0 or 1."

    N = data['X'].shape[0]
    D = data['X'].shape[1]

    # if Z does not exist, initialize
    if not(hidden.has_key('Z')):
        hidden['Z'] = 1.0*(np.random.rand(N,2) > 0.8)
        if params['bias'] == 1:
            hidden['Z'] = np.concatenate(np.ones((N,1)), hidden['Z'])

    # replace missings
    data['X'][np.isnan(data['X'])] = params['missing']

    # change labels for categorical and ordinal vars such that > 0
    V_offset = np.zeros(D)
    for d in xrange(D):
        if (data['C']=='c' or data['C']=='o'):
            mask = not(data['X'][:,d] == params['missing'])
            V_offset[d] = np.min( data['X'][mask,d] )
            data['X'][mask,d] = data['X'][mask,d] - V_offset[d] + 1

    # eventually, apply external transform
    for r in xrange(data['X'].shape[1]):
        if not(params['t'][r] == None): # there is an external transform
            data['X'][:,r] = params['t_1'][r](data['X'][:,r])
            data['C'][r] = params['ext_dataType'][r]

    # prepare input data for C++ inference routine
    Fin = np.ones(data['X'].shape[1]) # choose internal transform function (for positive)
    Xin = np.ascontiguousarray( data['X'].transpose() ) # specify way to store matrices to be
    Zin = np.ascontiguousarray( hidden['Z'].transpose() ) # compatible with C code
    tic = timeI.time()

    # RUN C++ routine
    (Z_out,B_out,Theta_out,mu_out,w_out,s2Y_out) = \
            GLFMlib.infer(Xin, data['C'], Zin, Fin, params['bias'], params['s2u'],\
            params['s2B'], params['alpha'], params['Niter'],\
            params['maxK'], params['missing'], params['verbose'])
    hidden['time'] = timeI.time() - tic
    if params['verbose']:
        print '\n\tElapsed time: %.2f seconds.\n' % hidden['time']

    # wrap output values inside hidden
    hidden['Z'] = Z_out.transpose()
    hidden['B'] = B_out
    hidden['theta'] = Theta_out
    hidden['mu'] = mu_out
    hidden['w'] = w_out
    hidden['s2Y'] = s2Y_out

    hidden['R'] = np.ones(D)
    for d in xrange(D):
        if (data['C']== 'c' or data['C'] == 'o'):
            hidden['R'][d] = np.unique( data['X']\
                    [not(data['X'][:,d] == params['missing']),d] ).shape[0]
    return hidden

def complete(data, hidden=dict(), params=dict()):
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
    # complete dictionary params
    params = init_default_params(data, params) # complete unspecified fields

    # check input syntax
    assert(type(data) is dict), "input 'data' should be a dictionary."
    assert(type(hidden) is dict), "input 'hidden' should be a dictionary."
    assert(type(params) is dict), "input 'params' should be a dictionary."
    assert(data.has_key('X')), "dictionary data does not have any matrix X defined."
    assert(data.has_key('C')), "dictionary data does not have any datatype array C defined."
    assert(params['bias'] <= 1), "bias parameter misspecified: should be either 0 or 1."

    pdb.set_trace()
    if sum( sum( (np.isnan(data['X'])) | (data['X']==params['missing']) )) == 0:
        print "The input matrix X has no missing values to complete."
        Xcompl = []
        return (Xcompl,hidden)

    # Run Inference
    hidden = infer(data,hidden,params)

    # Just in case there is any nan (also considered as missing)
    data['X'][np.isnan(data['X'])] = params['missing']

    [idxs_d, idxs_n] = (Xmiss == missing).nonzero()

    Xcompl=np.copy(Xmiss)
    for ii in xrange(len(idxs_n)): # for each missing
        if Xmiss[idxs_d[ii],idxs_n[ii]] == missing: # will always be the case
            Xcompl[xx_miss[i],yy_miss[i]] = computeMAP( C, Z[xx_miss[i],:], hidden, params, yy_miss[i] )
    return Xcompl

def computeMAP(C, Zp, hidden, params=dict(), idxsD=[]):
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

    if len(Zp.shape) == 1: # Zp is just 1 vector
        P = 1
        K2 = Zp.shape[0]
    else:
        P = Zp.shape[0]
        K2 = Zp.shape[1]
    K = hidden['B'].shape[1]
    assert (K2 == K), "Incongruent sizes between Zp and hidden['B']: number of latent variables should not be different"

    X_map = np.zeros((P,len(idxsD))) # output
    for dd in xrange(len(idxsD)): # for each dimension
        d = idxsD[dd]
        if params.has_key('t'): # if external transformations have been defined
            if not(params['t'][d] == None): # there is an external transform for data type d
                data['C'][d] = params['ext_dataType'][d] # set new type of data

        if not(data['C'][d] == 'c'):
            aux = np.inner(Zp, np.squeeze(hidden['B'][d,:]))

        if data['C'][d] == 'g':
            X_map[:,d] = mf.f_g( aux, hidden['mu'][d], hidden['w'][d] )
        elif data['C'][d] == 'p':
            pdb.set_trace()
            X_map[:,d] = f_p( aux, hidden['mu'][d], hidden['w'][d] )
        elif data['C'][d] == 'n':
            pdb.set_trace()
            X_map[:,d] = f_n( aux, hidden['mu'][d], hidden['w'][d] )
        elif data['C'][d] == 'c':
            pdb.set_trace()
            X_map[:,d] = f_c( np.inner(Zp, np.squeeze(hidden['B'][d,:,range(hidden['R'][d])])) )
        elif data['C'][d] == 'o':
            pdb.set_trace()
            X_map[:,d] = f_o( aux, hidden['theta'][d,range(hidden['R'][d]-1)] )
        else:
            raise ValueError('Unknown data type')
        if (sum(np.isnan(X_map[:,d])) > 0):
            raise ValueError('Some values are nan!')
        if params.has_key('t'):
            if not(params['t'][d] == None): # there is an external transform for data type d
                X_map[:,d] = params['t'][d]( X_map[:,d] )
    return X_map

def computePDF(data, Zp, hidden, params, d):
    """
    Function to compute probability density function for dimension d
    """
    data['X'][np.isnan(data['X'][:,d]),d] = params['missing_val']

    # compute x-domain [mm MM] to compute pdf
    mm = np.min(data['X'][not(data['X'][:,d] == params['missing']), d]) # min value
    MM = np.max(data['X'][not(data['X'][:,d] == params['missing']), d]) # max value

    if (not(params['t'][d]) == None): # if there is an external transformation
        data['C'][d] = params['ext_dataType'][d]
        mm = params['t_1'][d](mm)
        MM = params['t_1'][d](MM)

    if len(Zp.shape) == 1: # Zp is just 1 vector
        P = 1
        K2 = Zp.shape[0]
    else:
        P = Zp.shape[0]
        K2 = Zp.shape[1]
    K = hidden['B'].shape[1]
    assert (K2 == K), "Incongruent sizes between Zp and hidden['B']: number of latent variables should not be different"

    if (data['C'][d] == 'g') or (data['C'][d] == 'p'):
        if not(params.has_key('numS')):
            params['numS'] = 100
        xd = np.linspace(mm, MM, num=params['numS'])
    elif (data['C'][d] == 'n'):
        xd = range(mm,MM+1)
        params['numS'] = len(xd)
    else:
        xd = np.unique(data['X'][not(data['X'][:,d] == params['missing']), d])
        params['numS'] = len(xd) # number of labels for categories or ordinal data
    pdf = np.zeros((P,params['numS']))

    for p in xrange(P):
        if data['C'][d] == 'g':
            pdf[p,:] = mf.pdf_g(xd,Zp[p,:], np.squeeze(hidden['B'][d,:]), hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d], params)
        elif data['C'][d] == 'p':
            pdf[p,:] = mf.pdf_p(xd,Zp[p,:], np.squeeze(hidden['B'][d,:]), hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d], params)
        elif data['C'][d] == 'n':
            pdf[p,:] = mf.pdf_n(xd,Zp[p,:], np.squeeze(hidden['B'][d,:]), hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d], params)
        elif data['C'][d] == 'c':
            pdf[p,:] = mf.pdf_c(Zp[p,:], np.squeeze(hidden['B'][d,:,range(hidden['R'][d])]), hidden['s2Y'][d])
        elif data['C'][d] == 'o':
            pdf[p,:] = mf.pdf_o(Zp[p,:], squeeze(hidden['B'][d,:]), hidden['theta'][d,range(hidden['R'][d]-1)], hidden['s2Y'][d])
        else:
            raise ValueError('Unknown data type')
        assert (np.sum(np.isnan(pdf)) == 0), "Some values are nan!"

    if params.has_key('t'):
        if not(params['t'] == None): # we have used a special transform beforehand
            xd = params['t'][d](xd) # if there was an external transformation, transform pdf
            pdf = pdf * np.abs( params['dt_1'][d](xd) )
    return (xd,pdf)

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
    Xd = data['X'][d,:]
    data['C'] = data['C'][d]
    if k<0 or k>Kest:
        print('Error: k index should be bigger than o and smaller than Kest')
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    if data['C'] == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
        plt.plot(xx,pdf)

    elif data['C'] == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w));
        plt.plot(xx,pdf)

    elif data['C'] == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
        plt.stem(xx,pdf)

    elif data['C'] == 'c':
        R = len( np.unique(Xd[mask]) )
        Zn = np.zeros(Kest)
        Zn[k] = 1
        Bdv = np.squeeze(B[d,:,:]) # TODO: Review that size = [K*R]
        pdf = mf.pdf_cat(Zn,Bdv,s2u,R)
        bar_width = 0.35
        index = np.arange(len(pdf))
        plt.bar(index,pdf,bar_width)
        plt.xticks(index + bar_width / 2, catlabel, rotation='vertical')

    elif data['C'] == 'o':
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
    Xd = data['X'][d,:]
    data['C'] = data['C'][d]
    # only consider values in dimension d which are not missing
    if np.isnan(missing):
        mask = np.where(not(np.isnan(Xd)))[0]
    else:
        mask = np.where(Xd != missing)[0]
    (numPatterns,Kest) = Zp.shape
    if data['C'] == 'g':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_real(xx, Zn,Bdv,s2Y,s2u)
            plt.plot(xx,pdf,label=str(Zn))

    elif data['C'] == 'p':
        numP = 100 # number of points to plot
        xx = np.linspace( min(Xd[mask]), max(Xd[mask]), numP )
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_pos(xx,Zn,Bdv,w,s2Y,s2u,lambda x,w: mf.fpos_1(x,w), \
                lambda x,w: mf.dfpos_1(x, w))
            plt.plot(xx,pdf,colors[p],label=str(Zn))

    elif data['C'] == 'n':
        xx = np.arange( min(Xd[mask]), max(Xd[mask])+1)
        Bdv = B[d,:,0]
        w = 2.0 / max(Xd[mask]) # TODO: put function handler
        for p in xrange(numPatterns):
            Zn = np.squeeze(Zp[p,:]) # TODO: Verify dimensions
            pdf = mf.pdf_count(xx,Zn,Bdv,w,s2Y, lambda x,w: mf.fpos_1(x,w))
            plt.stem(xx,pdf,colors[p], label=str(Zn))

    elif data['C'] == 'c':
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

    elif data['C'] == 'o':
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


def init_default_params(data, params):
    """
    Initialize of complete dictionary params
    Input:
        data: dict with database
        params: dict to complete with default values
    Output:
        same data structure: params
    """
    # s2u=0.001
    D = data['X'].shape[1]
    if not(params.has_key('missing')):
        params['missing'] = -1
    if not(params.has_key('alpha')):
        params['alpha'] = 1
    if not(params.has_key('bias')):
        params['bias'] = 0
    if not(params.has_key('s2u')):
        params['s2u'] = 0.01
    if not(params.has_key('s2B')):
        params['s2B'] = 1
    if not(params.has_key('Niter')):
        params['Niter'] = 1000
    if not(params.has_key('maxK')):
        params['maxK'] = D
    if not(params.has_key('verbose')):
        params['verbose'] = 1
    if not(params.has_key('numS')):
        params['numS'] = 1

    # parameters for optional external transformation
    if not(params.has_key('t')):
        params['t'] = [None] * D
    if not(params.has_key('t_1')):
        params['t_1'] = [None] * D
    if not(params.has_key('dt_1')):
        params['dt_1'] = [None] * D
    return params
