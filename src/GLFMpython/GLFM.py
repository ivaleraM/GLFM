import numpy as np
import random

import mapping_functions as mf
import matplotlib.pyplot as plt

import time
import os
import sys
root = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])
sys.path.append(os.path.join(root, 'Ccode/wrapper_python/'))

import GLFMlib # python wrapper library in order to run C++ inference routine
import mapping_functions as mf

import copy

def infer(data,hidden=dict(), params=dict()):
    """
    Python wrapper to launch inference routine for the GLFM model.
    Inputs:
        data: dictionary containing all input data structures
            X: input data N*D where:
                N = number of observations
                D = number of dimensions
            C: string array indicating types of data ('g': real,'p': positive real,
                'c': categorical; 'o': ordinal; 'n': count data)
        hidden (optional): dictionary containing latent variables
            Z: initial feature activation matrix: N*K
                K = number of latent dimensions
        params (optional): dictionary containing eventual simulation parameters
            bias: indicator of whether to include or not a bias
            s2u: internal auxiliary noise
            s2B: noise variance for prior over elements of matrix B
            alpha: concentration parameter of the Indian Buffet Process
            Niter: number of simulations
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
    # complete dictionary params with default values
    params = init_default_params(data, params) # complete unspecified fields

    # check input syntax
    assert(type(data) is dict), "input 'data' should be a dictionary."
    assert(type(hidden) is dict), "input 'hidden' should be a dictionary."
    assert(type(params) is dict), "input 'params' should be a dictionary."
    assert(data.has_key('X')), "dictionary data does not have any matrix X defined."
    assert(data.has_key('C')), "dictionary data does not have any datatype array C defined."
    assert(params['bias'] <= 1), "bias parameter misspecified: should be either 0 or 1."

    N = data['X'].shape[0] # number of observations
    D = data['X'].shape[1] # number of dimensions

    # if Z does not exist, initialize
    if not(hidden.has_key('Z')):
        hidden['Z'] = 1.0*(np.random.rand(N,2) > 0.8)
        if params['bias'] == 1: # add bias if requested
            hidden['Z'] = np.concatenate((np.ones((N,1)), hidden['Z']),axis=1)

    tmp_data = copy.deepcopy(data) #np.copy(data).tolist()

    # replace nan by missing values
    tmp_data['X'][np.isnan(tmp_data['X'])] = params['missing']
    # # dealing with missing data: replace np.nan by -1
    # (xx,yy) = np.where(np.isnan(X)) # find positions where X is nan (i.e. missing data)
    # for r in xrange(len(xx)):
    #     X[xx[r],yy[r]] = -1

    # change labels for categorical and ordinal vars such that categories
    # start counting at 1 and all of them are bigger than 0
    V_offset = np.zeros(D)
    for d in xrange(D):
        if (tmp_data['C'][d]=='c' or tmp_data['C'][d]=='o'):
            mask = tmp_data['X'][:,d] != params['missing']
            uniqueVal = np.unique(tmp_data['X'][mask,d])
            Xaux = np.zeros(N)
            for i in xrange(len(uniqueVal)):
                Xaux[tmp_data['X'][:,d] == uniqueVal[i]] = i+1
            Xaux[map(lambda x: not x, mask)] = params['missing']
            tmp_data['X'][:,d] = Xaux
            #V_offset[d] = np.min( tmp_data['X'][mask,d] )
            #tmp_data['X'][mask,d] = tmp_data['X'][mask,d] - V_offset[d] + 1

    # eventually, apply external transform specified by the user
    for r in xrange(tmp_data['X'].shape[1]):
        if not(params['t'][r] == None): # there is an external transform
            tmp_data['X'][:,r] = params['t_1'][r](tmp_data['X'][:,r])
            tmp_data['C'] = tmp_data['C'][:r] + params['ext_dataType'][r] + tmp_data['C'][(r+1):]

    # prepare input data for C++ inference routine
    Fin = np.ones(tmp_data['X'].shape[1]) # choose internal transform function (for positive)
    Xin = np.ascontiguousarray( tmp_data['X'].transpose() ) # specify way to store matrices to be
    Zin = np.ascontiguousarray( hidden['Z'].transpose() ) # compatible with C code
    tinit = time.time() # start counting time

    # RUN C++ routine
    (Z_out,B_out,Theta_out,mu_out,w_out,s2Y_out) = \
            GLFMlib.infer(Xin, tmp_data['C'], Zin, Fin, params['bias'], params['s2u'],\
            params['s2B'], params['alpha'], params['Niter'],\
            params['maxK'], params['missing'], params['verbose'])

    tlast = time.time()

    hidden['time'] = tlast - tinit
    if params['verbose']:
        print '\n\tElapsed time: %.2f seconds.\n' % hidden['time']

    # wrap output values inside hidden
    hidden['Z'] = Z_out.transpose()
    hidden['B'] = B_out
    hidden['theta'] = Theta_out
    hidden['mu'] = mu_out
    hidden['w'] = w_out
    hidden['s2Y'] = s2Y_out

    hidden['R'] = [None] * D
    for d in xrange(D):
        if (data['C'][d] == 'c' or data['C'][d] == 'o'):
            hidden['R'][d] = np.unique( data['X']\
                    [data['X'][:,d] != params['missing'],d] )
    return hidden

def complete(data, hidden=dict(), params=dict()):
    """
    Inputs:
        data: dictionary containing all input data structures
            X: input data N*D where:
                N = number of observations
                D = number of dimensions
            C: string array indicating types of data ('g': real,'p': positive real,
                'c': categorical; 'o': ordinal; 'n': count data)
        hidden (optional): dictionary containing latent variables
            Z: initial feature activation matrix: N*K
                K = number of latent dimensions
        params (optional): dictionary containing eventual simulation parameters
            bias: indicator of whether to include or not a bias
            s2u: internal auxiliary noise
            s2B: noise variance for prior over elements of matrix B
            alpha: concentration parameter of the Indian Buffet Process
            Niter: number of simulations
            maxK: maximum number of features for memory allocation
            missing: value for missings (should be an integer, not nan)
            verbose: indicator to print more information
    Output:
        Xcompl : same numpy array as input X whose missing values have been
                 inferred and completed by the algorithm.
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

    if sum( sum( (np.isnan(data['X'])) | (data['X']==params['missing']) )) == 0:
        print "The input matrix X has no missing values to complete."
        Xcompl = []
        return (Xcompl,hidden)

    # Run Inference
    hidden = infer(data,hidden,params)

    tmp_data = copy.deepcopy(data) # np.copy(data).tolist()
    # Just in case there is any nan (also considered as missing)
    tmp_data['X'][np.isnan(tmp_data['X'])] = params['missing']

    [xx_miss, yy_miss] = (tmp_data['X'] == params['missing']).nonzero()

    Xcompl=np.copy(tmp_data['X'])
    for ii in xrange(len(xx_miss)): # for each missing
        if tmp_data['X'][xx_miss[ii],yy_miss[ii]] == params['missing']: # will always be the case
            Xcompl[xx_miss[ii],yy_miss[ii]] = computeMAP( tmp_data['C'], hidden['Z'][xx_miss[ii],:], hidden, params, [ yy_miss[ii] ] )
    return (Xcompl,hidden)

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
          - s2Y: 1*D inferred noise variance for each dimension of pseudo-observations
          - theta: D*maxR matrix of auxiliary vars (for ordinal variables)
    ----------------(optional) ------------------
          - idxsD: dimensions to infer

    Outputs:
      X_map: P*Di matrix with MAP estimate where Di = length(idxsD)
    """
    if (len(idxsD) == 0): # no dimension specified, infer all dimensions
        idxsD = range(hidden['B'].shape[0])

    if len(Zp.shape) == 1: # Zp is just 1 vector
        P = 1
        K2 = Zp.shape[0]
    else:
        P = Zp.shape[0]
        K2 = Zp.shape[1]
    K = hidden['B'].shape[1] # number of latent features
    assert (K2 == K), "Incongruent sizes between Zp and hidden['B']: number of latent variables should not be different"

    X_map = np.zeros((P,len(idxsD))) # output matrix
    for dd in xrange(len(idxsD)): # for each dimension
        d = idxsD[dd]
        if params.has_key('t'): # if external transformations have been defined
            if not(params['t'][d] == None): # there is an external transform for data type d
                C = C[:d] + params['ext_dataType'][d] + C[(d+1):]

        if not(C[d] == 'c'):
            aux = np.inner(Zp, hidden['B'][d,:,0])

        if C[d] == 'g':
            X_map[:,dd] = mf.f_g( aux, hidden['mu'][d], hidden['w'][d] )
        elif C[d] == 'p':
            X_map[:,dd] = mf.f_p( aux, hidden['mu'][d], hidden['w'][d] )
        elif C[d] == 'n':
            X_map[:,dd] = mf.f_n( aux, hidden['mu'][d], hidden['w'][d] )
        elif C[d] == 'c':
            X_map[:,dd] = mf.f_c( np.inner(Zp, hidden['B'][d,:,\
                    range(int(hidden['R'][d].shape[0])) ]) )
        elif C[d] == 'o':
            X_map[:,dd] = mf.f_o( aux, hidden['theta'][d,range(int(hidden['R'][d].shape[0]-1))] )
        else:
            raise ValueError('Unknown data type')
        if (sum(np.isnan(X_map[:,dd])) > 0):
            raise ValueError('Some values are nan!')
        if params.has_key('t'):
            if not(params['t'][d] == None): # there is an external transform for data type d
                X_map[:,dd] = params['t'][d]( X_map[:,dd] )
    return X_map

def computePDF(data, Zp, hidden, params, d):
    """
    Function to compute probability density function for dimension d
    """
    print "dim=%d\n" % d
    tmp_data = copy.deepcopy(data) #np.copy(data).tolist()
    tmp_data['X'][np.isnan(tmp_data['X'][:,d]),d] = params['missing']

    # compute x-domain [mm MM] to compute pdf
    mm = np.min(tmp_data['X'][tmp_data['X'][:,d] != params['missing'], d]) # min value
    MM = np.max(tmp_data['X'][tmp_data['X'][:,d] != params['missing'], d]) # max value

    if (params['t'][d] != None): # if there is an external transformation
        tmp_data['C'] = tmp_data['C'][:d] + params['ext_dataType'][d] + tmp_data['C'][d+1:]
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

    if (tmp_data['C'][d] == 'g') or (tmp_data['C'][d] == 'p'):
        numS = 100
        xd = np.linspace(mm, MM, num=numS)
    elif (tmp_data['C'][d] == 'n'):
        xd = np.array( range(int(mm),int(MM)+1) )
        numS = len(xd)
    else:
        xd = np.unique(tmp_data['X'][tmp_data['X'][:,d] != params['missing'], d])
        numS = len(xd) # number of labels for categories or ordinal data
    pdf = np.zeros((P,numS))

    for p in xrange(P):
        if tmp_data['C'][d] == 'g':
            pdf[p,:] = mf.pdf_g(xd,Zp[p,:], hidden['B'][d,:,0], hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d], params['s2u'])
        elif tmp_data['C'][d] == 'p':
            pdf[p,:] = mf.pdf_p(xd,Zp[p,:], hidden['B'][d,:,0], hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d], params['s2u'])
        elif tmp_data['C'][d] == 'n':
            pdf[p,:] = mf.pdf_n(xd,Zp[p,:], hidden['B'][d,:,0], hidden['mu'][d], hidden['w'][d], hidden['s2Y'][d])
        elif tmp_data['C'][d] == 'c':
            pdf[p,:] = mf.pdf_c(Zp[p,:], hidden['B'][d,:,range(int(hidden['R'][d].shape[0]))], hidden['s2Y'][d])
        elif tmp_data['C'][d] == 'o':
            pdf[p,:] = mf.pdf_o(Zp[p,:], hidden['B'][d,:,0], hidden['theta'][d,range(int(hidden['R'][d].shape[0]-1))], hidden['s2Y'][d])
        else:
            raise ValueError('Unknown data type')
        assert (np.sum(np.isnan(pdf)) == 0), "Some values are nan!"

    if params.has_key('t'):
        if (params['t'][d] != None): # we have used a special transform beforehand
            xd = params['t'][d](xd) # if there was an external transformation, transform pdf
            pdf = pdf * np.abs( params['dt_1'][d](xd) )
    return (xd,pdf)

def get_feature_patterns_sorted(Z):
    """
    Function to compute list of activation patterns. Returns sorted list
    Input:
        Z: N*K binary matrix
    Outputs:
        patterns: numP*K: list of patterns
        C: assignment vector of length N*1 with pattern id for each observation
        L: numP*1 vector with num. of observations per pattern
    """
    N = Z.shape[0]
    C = np.zeros(N)

    patterns = np.vstack({tuple(row) for row in Z})
    numP = patterns.shape[0]
    L = np.zeros(numP)
    for r in xrange(numP): # for each pattern
        pat = patterns[r,:]
        mask = np.sum(np.tile(pat,(N,1)) == Z, axis=1) == Z.shape[1]
        C[mask] = r
        L[r] = sum(mask)
        #print '%d. %s: %d' % (r, str(patterns[r,:]), L[r])

    # sort arrays
    idxs = L.argsort()
    idxs = np.flipud(idxs)
    L = L[idxs]
    C = C[idxs]
    patterns = patterns[idxs,:]

    print '\n'
    for r in xrange(numP): # for each pattern
        print '%d. %s: %d' % (r, str(patterns[r,:]), L[r])

    return (patterns,C,L)

def plotPatterns(data, hidden, params, patterns, colors=[], styles=[],\
        leg=[], idxD=[]):
    """
    Function to plot the inferred distribution and empirical histogram for
    each dimension of the observations.
    Inputs:
        data: data structure
        hidden: structure of latent variables
        params: structure of simulation parameters and hyperparameters
        patterns: numP*K list of patterns to plot
        ------ (optional) ------
        colors: list of colors to plot
        styles: list of styles for each line (for plot, not bar)
        leg: legend to use (by default, use patterns as legend)
        idxD: array of dimensions to plot
    Outputs:
        void
    """

    # initialize optional parameters
    if len(leg) == 0: # use patterns as legend
        leg = [str(x) for x in patterns ]
    if len(idxD) == 0:
        idxD = range(data['X'].shape[1])

    if (patterns.shape[1] != hidden['B'].shape[1]):
        raise ValueError('Error: Sizes of patterns and B are inconsistent')
    # if (len(leg) != patterns.shape[0]):
    #     raise ValueError('Error: Sizes of leg and patterns are inconsistent')

    (D,Kest,maxR) = hidden['B'].shape
    (numPatterns,Kest) = patterns.shape

    # replace nan by missing values
    tmp_data = copy.deepcopy(data) #np.copy(data).tolist()
    tmp_data['X'][np.isnan(tmp_data['X'])] = params['missing']

    colors = [[0.8784, 0.8784, 0.8784], 'r','b','g','m','k', \
            [0.9290, 0.6940, 0.1250], [0.4660, 0.6740, 0.1880]] # default colors

    # legend for empirical histogram
    leg = ['Empirical'] + leg

    for dd in xrange(len(idxD)): # for each required dimension
        d = idxD[dd]
        plt.figure(d)                 # create new figure
        #hold off # TODO
        #plt.xlabel(data['ylabel'][d]) # add x legend
        # plot empirical if 'g' | 'p' | 'n'
        if (data['C'][d] == 'g' or data['C'][d] == 'p' or data['C'][d] == 'n'):
            mask = tmp_data['X'][:,d]!=params['missing']
            #plt.hist(data['X'][ mask ,d], bins='auto')
            results, edges = np.histogram(tmp_data['X'][ mask ,d], normed=True, bins=100)
            binWidth = edges[1] - edges[0]
            plt.bar(edges[:-1], results, binWidth, color=[0.7529, 0.7529, 0.7529])
            plt.show()

          #  ax = plt.subplot(111)
          #  w = 0.9 / float(Kest)
          #  ax.bar(x-w, y,width=w,color='b',align='center')
          #  ax.bar(x, z,width=w,color='g',align='center')
          #  ax.bar(x+w, k,width=w,color='r',align='center')
          #  ax.xaxis_date()
          #  ax.autoscale(tight=True)

            #[h xx] = hist(data.X(:,d),100);
            #h = h ./ sum(h * (xx(2) - xx(1)));
            #bar(xx, h);
            #set(get(gca,'child'),'FaceColor',[0.8784 0.8784 0.8784], ...
            #    'EdgeColor',[0.7529 0.7529 0.7529]);
            # #hold on;

        (xd,pdf) = computePDF(tmp_data, patterns, hidden, params, d)
        if (tmp_data['C'][d] == 'c') or (tmp_data['C'][d] == 'o'):
            mask = tmp_data['X'][:,d] != params['missing']
            (tmp,bla) =  np.histogram(tmp_data['X'][:,d], \
                   np.unique(tmp_data['X'][:,d]).tolist() + \
                   [np.unique(tmp_data['X'][:,d])[-1] + 1],\
                   density=True)
            # tmp = hist(data.X(mask,d), unique(data.X(mask,d)));
            bar_width = 0.8/(numPatterns+1)
            plt.bar(xd,tmp,width=bar_width, color=colors[0], label=leg[0]) # plot empirical
            for p in xrange(numPatterns):
                plt.bar(xd+(p+1)*bar_width,pdf[p,:],width=bar_width,\
                        color=colors[np.mod(p+1,len(colors))], label=leg[p+1])
            plt.xticks(xd + 0.4, data['cat_labels'][d]) #, rotation='vertical')

            #tmp = tmp / sum(tmp);
            #h = bar([tmp' pdf']);
            #h(1).FaceColor = [0.8784 0.8784 0.8784];

        elif (tmp_data['C'][d] == 'n'):
            for p in xrange(numPatterns):
                markerline, stemlines, baseline = plt.stem(xd,pdf[p,:], label=leg[p+1])
                plt.setp(stemlines, 'color', colors[np.mod(p+1,len(colors))]) # plt.getp(markerline,'color'))
                plt.setp(markerline, 'markerfacecolor', colors[np.mod(p+1,len(colors))]) # plt.getp(markerline,'color'))

        else:
            #inte = np.zeros(pdf.shape[0])
            #for p in xrange(numPatterns):
            #    for r in xrange(len(xd)):
            #        inte[p] = inte[p] + (xd[p+1] - xd[p])*pdf[p,r]
            for p in xrange(numPatterns):
            #    print "int = %f\n" % sum((xd[1]-xd[0])*pdf[p,:])
                plt.plot(xd,pdf[p,:], color=colors[np.mod(p+1,len(colors))], label=leg[p+1])

        #if len(colors) > 0:
        #   # set colors in plot
        #
        #if len(styles) > 0:
        #   # set styles in plot

        plt.legend()
        plt.ion()
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

    # parameters for optional external transformation
    if not(params.has_key('t')):
        params['t'] = [None] * D
    if not(params.has_key('t_1')):
        params['t_1'] = [None] * D
    if not(params.has_key('dt_1')):
        params['dt_1'] = [None] * D
    return params
