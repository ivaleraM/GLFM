# -------------------------------------------------------------------------
# Module with all available mapping functions necessary for the Generalized
# Latent Feature Model (GLFM)
# -------------------------------------------------------------------------
import numpy as np
from scipy.stats import norm

import pdb

# --------------------
# Mapping Functions
# --------------------

def f_g(y, mu, w):
    # Mapping function for real-valued data
    #  Y -> X (from pseudo-obversations to data)
    assert not(w == 0), 'scaling factor should never be 0'
    x = (y*1.0)/w + mu;
    return x

def f_p(y, mu, w):
    # transformation function for positive data
    # Y -> X (from pseudo-obversations to data)
    assert not(w == 0), 'scaling factor should never be 0'
    x = np.log( np.exp(y) + 1 )*1.0/w + mu
    return x

def f_c(y):
    # transformation function for categorical data
    # input argument y: [N*R]
    # output: x [N*1]
    assert (len(y.shape) > 1), 'there is only one category, this dimension does not make sense'
    x = np.zeros(y.shape[0])
    for i in xrange(y.shape[0]):
        val = max(y[i,:])
        x[i] = np.where(y[i,:] == val)[0][0] + 1
    # [a,x] = max(y,[],2);
    return x

def f_n(y,mu,w):
    # transformation function for count data
    # Y -> X (from pseudo-obversations to data)
    # Inputs:
    #   Pseudo-observations y: [N*D]
    assert not(w == 0), 'scaling factor should never be 0'
    x = np.floor( f_p(y,mu,w) )
    return x

def f_o(y, theta):
    # Mapping function for ordinal data
    # Inputs:
    #       y: [1*R] Pseudo-observations
    #   theta: [1*(R-1)] Thresholds that divide the real line into R regions
    x = np.zeros(y.shape[0]) # column vector
    for j in xrange(len(theta)):
        if (j == 0):
            mask = (y <= theta[0])
        else:
            mask = (y > theta[j-1]) * (y <= theta[j])
        x[mask] = j + 1
    x[x == 0] = len(theta) # last ordinal category
    #for r in xrange(len(theta)): #= 1:length(theta)
    #    val = theta[r]
    #    if (y < val):
    #        x = r
    #        break
    #    if r == size(theta,3):
    #        x = r+1
    return x

def f_g_1(x, mu, w):
    # transformation function for real-valued data
    # X -> Y (from data to pseudo-obversations)
    assert not(w == 0), 'scaling factor should never be 0'
    y = w * (x - mu)
    return y

def f_n_1(x, mu, w):
    # transformation function for positive data
    # X -> Y (from data to pseudo-obversations)
    y = f_p_1(x, mu, w)
    return y

def f_p_1(x, mu, w):
    # transformation function for positive data
    # X -> Y (from data to pseudo-obversations)
    assert not(w == 0), 'scaling factor should never be 0'
    try:
        y = np.log( np.exp(w*(x-mu) - 1) )
        #y = w*(x-mu) - 1
    except:
        print('Waiting at debugging point')
        pdb.set_trace()
    return y

def df_p_1(x, mu, w):
    # derivative of transformation function for positive data
    # X -> Y (from data to pseudo-obversations)
    assert not(w == 0), 'scaling factor should never be 0'
    try:
        y = ( w * np.exp(w*(x-mu)) ) / ( np.exp(w*(x-mu) - 1) )
    except:
        print('Waiting at debugging point')
        pdb.set_trace()
    return y

# ------------------------------------------
# Functions to compute pdf values
# ------------------------------------------

def pdf_g(X,Zp, Bd, mu , w, s2y, s2u):
    """
    Probability Density Function for real variables
    Eq. (2) in the paper:
    "General Latent Feature Models for Heterogeneous Datasets"

    Inputs:
      X: univariate observations (for a particular dimension)
      Zp: 1*K array
      Bd: K*D array
     s2y: noise variance for pseudo-observations (scalar)
     s2u: auxiliary noise variance (scalar)
    """
    df_1 = lambda x: w*(x-mu)
    pdf = norm.pdf( df_1(X) , np.dot(Zp,Bd) , np.sqrt(s2y + s2u)) * w
    return pdf

def pdf_c(Zn,B,s2y,numMC_samples=100):
    """
    Function to compute pdf of a categorical variable. It returns the whole
    pdf, a probability vector of length R (number of categories)
    Eq. (4) in the paper:
    "General Latent Feature Models for Heterogeneous Datasets"

    Inputs:
      B: [R*K]  feature weight matrix (dictionary)
          where R: number of categories
      Zn: [K], binary vector of feature assignment,
          where K: number of latent features
      s2y: scalar, variance of auxiliary noise
    """

    R = B.shape[0]
    pdf = np.zeros(R)
    uV = np.sqrt(s2y) * np.random.randn(numMC_samples) # mean for u = 0
    for r in xrange(R):
        tmp = np.ones((1,numMC_samples))
        # we compute the expectation using Monte Carlo samples
        for j in xrange(R):
            if (j==r):
                continue
            tmp = tmp * norm.cdf(uV + np.kron( np.ones(numMC_samples), \
                    np.inner(Zn,B[r,:]-B[j,:]) ) )
            #aux = aux .* normcdf( u + Zp*(B(:,r) - B(:,j)) );
        pdf[r] = np.mean(tmp,1)
    pdf = pdf / sum(pdf)
    return pdf

def pdf_n(X,Zn,Bd,mu,w,s2y):
    """
    Likelihood function for count data
    Eq. (8) in the paper:
    "General Latent Feature Models for Heterogeneous Datasets"

    Inputs:
       X: count-data observation
      Bd: K*1  weight vector for a particular dimension d
      Zn: 1*K, binary vector of feature assignment,
          where K: number of latent features
     s2y: Gaussian noise variance  for u
    mu,w: Hyper-parameters linked to the transformation described in the paper
    """

    pdf = norm.cdf(f_n_1(X+1,mu,w), np.inner(Zn,Bd), np.sqrt(s2y)) \
        - norm.cdf(f_n_1(X,mu,w), np.inner(Zn,Bd), np.sqrt(s2y))
    return pdf

def pdf_o(Zn,Bd,theta,s2y):
    """
    Function to compute pdf of an ordinal variable. It returns the whole
    pdf, a probability vector of length R (number of categories)
    Eq. (6) in the paper:
    "General Latent Feature Models for Heterogeneous Datasets"

    Input parameters:
       Zn: [1*K], binary feature activation vector
       Bd: [K*1], feature weight vector (dictionary) for dimension d
    theta: [1*(R-1)]
      s2y: scalar, variance of pseudo-observations
    theta: 1xR auxiliary thresholds, see description in the paper
    """

    R = len(theta)+1 # number of categories
    pdf = np.zeros(R)
    for r in xrange(R):
        pdf[r] = pdf_o_single(r,Zn,Bd,theta,s2y)
    return pdf

def pdf_o_single(r,Zn,Bd,theta,s2y):
    """
    Function to compute p(X=r|Z,B) for ordinal value r
    """
    R = len(theta)+1 # number of categories
    assert (r>=0) and (r<= (len(theta)+1))
    if (r==0):
        a = norm.cdf(theta[0],np.inner(Zn,Bd),np.sqrt(s2y))
        b = 0
    elif (r==(R-1)):
        a = 1
        b = norm.cdf(theta[r-1],np.inner(Zn,Bd),np.sqrt(s2y))
    else:
        a = norm.cdf(theta[r],np.inner(Zn,Bd),np.sqrt(s2y))
        b = norm.cdf(theta[r-1],np.inner(Zn,Bd),np.sqrt(s2y))
    return (a-b)


def pdf_p(X,Zn,Bd,mu, w, s2y,s2u):
    """
    Probability density function for positive real variables
    Eq. (2) in the paper:
    "General Latent Feature Models for Heterogeneous Datasets"

    Inputs:
      X: values in which to evaluate pdf
     Zn: [1*K] vector, binary vector of feature assignment
     Bd: [K*1] weight vector for dimension d
     mu,w: scalar, hyper-parameters for transformation described in the paper
          In particular, w: normalization weight
          (this is computed by GLFM infer function and passed to this function)
    s2y: scalar, variance of pseudo-observations
    s2u: scalar, variance of auxiliary noise
    """

    pdf = 1.0/np.sqrt(2*np.pi*(s2u+s2y)) * np.exp( \
            -1.0/(2.0*(s2y+s2u)) * (f_p_1(X,mu,w) - np.inner(Zn,Bd))**2)\
            * np.abs( df_p_1(X,mu, w) )
    return pdf

## --------------------
## Auxiliary Functions
## --------------------
#
#def exp_count(Zn,Bd, s2y, fpos_1_handler, w, maxX):
#    # function to compute expectation of random variable
#    # Zn: 1*K nparray
#    # Bd: K*1 nparray
# #   s2y = params{1}
# #   fpos_1_handler = params{2}
# #   w = params{3} # function to compute normalization weights w
# #   maxX = params{4}
#    xin = range(1,maxX+1)
#    func = lambda x: pdf_count(x,Zn,Bd,w, s2y,fpos_1_handler)
#    expect = np.inner(xin,func(xin))
#    return expect
#
#def rnd_pos(Zn,Bd,numSamples,s2y,s2u,fpos_1_handler,dfpos_1_handler,w,maxX):
#    # function to get random samples from the distribution
#
#  #  s2y = params{1}
#  #  s2u = params{2}
#  #  fpos_1_handler = params{3}
#  #  dfpos_1_handler = params{4}
#  #  w = params{5} # function to compute normalization weights w
#  #  maxX = params{6}
#
#    func = lambda x: np.log( pdf_pos(x,Zn,Bd,w, s2y ,s2u, fpos_1_handler, dfpos_1_handler) )
#    a = 10**-6
#    b = maxX
#    domain = [a,b+5]
#
#    # FIND ARS IN PYTHON
#    samples = ars(func, a, b, domain, numSamples, [])
#    samples = np.exp(samples)
#    return samples
#
#def rnd_real(Zn,Bd,numSamples, s2y, s2u):
#    # function to get random samples from the distribution
#    x = np.sqrt(s2u+s2y) * np.random.randn(numSamples) + np.inner(Zn,Bd)
#    return x
