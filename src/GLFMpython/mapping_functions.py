# -------------------------------------------------------------------------
# Module with all available mapping functions necessary for the Generalized
# Latent Feature Model (GLFM)
# -------------------------------------------------------------------------
import numpy as np
from scipy.stats import norm

import pdb

# --------------------
# Auxiliary Functions
# --------------------

def exp_count(Zn,Bd, s2y, fpos_1_handler, w, maxX):
    # function to compute expectation of random variable
    # Zn: 1*K nparray
    # Bd: K*1 nparray
 #   s2y = params{1}
 #   fpos_1_handler = params{2}
 #   w = params{3} # function to compute normalization weights w
 #   maxX = params{4}
    xin = range(1,maxX+1)
    func = lambda x: pdf_count(x,Zn,Bd,w, s2y,fpos_1_handler)
    expect = np.inner(xin,func(xin))
    return expect

def rnd_pos(Zn,Bd,numSamples,s2y,s2u,fpos_1_handler,dfpos_1_handler,w,maxX):
    # function to get random samples from the distribution

  #  s2y = params{1}
  #  s2u = params{2}
  #  fpos_1_handler = params{3}
  #  dfpos_1_handler = params{4}
  #  w = params{5} # function to compute normalization weights w
  #  maxX = params{6}

    func = lambda x: np.log( pdf_pos(x,Zn,Bd,w, s2y ,s2u, fpos_1_handler, dfpos_1_handler) )
    a = 10**-6
    b = maxX
    domain = [a,b+5]

    # FIND ARS IN PYTHON
    samples = ars(func, a, b, domain, numSamples, [])
    samples = np.exp(samples)
    return samples

def rnd_real(Zn,Bd,numSamples, s2y, s2u):
    # function to get random samples from the distribution
    x = np.sqrt(s2u+s2y) * np.random.randn(numSamples) + np.inner(Zn,Bd)
    return x

# --------------------
# Mapping Functions
# --------------------

def freal(y, s2u):
    # Mapping function for real values
    # s2u: auxiliary noise
    dim = len(y.shape)
    if dim == 1:
        tmp = randn(y.shape[0])
    elif dim == 2:
        tmp = randn(y.shape[0],y.shape[1])
    else:
        print 'undefined dimensions for y'
    x = y + np.sqrt(s2u) * tmp
    return x

def fpos(y):
    # Inputs:
    #   Pseudo-observations y: [N*D] # Check dimensions
    x = np.log(np.exp(y)+1)
    return x

def fcat(y):
    # input argument y: [N*R]
    # output: x [1*N]
    x = np.max(y,1) # .reshape(-1,1)
    return x

def fcount(y):
    # Inputs:
    #   Pseudo-observations y: [N*D]
    x = np.floor( fpos(y) )
    return x

def ford(y, theta):
    # Mapping function for ordinal data
    # Inputs:
    #       y: [1*R] Pseudo-observations
    #   theta: [1*(R-1)] Thresholds that divide the real line into R regions
    for r in xrange(len(theta)): #= 1:length(theta)
        val = theta[r]
        if (y < val):
            x = r
            break
        if r == size(theta,3):
            x = r+1
    return x

def fpos_1(x):
    y = np.log( np.exp(x) - 1 )
    return y

def fpos_1_sq(x):
    y = np.sqrt(x)
    return y

def fpos_sq(y):
    # Inputs:
    #   y: pseudo-observations
    x =  y**2
    return x

def dfpos_1(x):
    y = 1.0 / ( 1 - np.exp(-x) )
    return y

def dfpos_1_xi(x, w):
    y = 2*x
    #y = -0.5*(w**0.5) * (x**(-0.5))
    return y

# ------------------------------------------
# Functions to compute pdf values
# ------------------------------------------

def pdf_real(X, Zn,Bd,s2y,s2u):
    # Probability Density Function for real variables
    # Inputs:
    #   X: observations (dimensions?)
    #   Zn: 1*K array
    #   Bd: 1*K array
    #  s2y: noise variance for pseudo-observations (scalar)
    #  s2u: auxiliary noise variance (scalar)
    pdf = norm.pdf(X, np.inner(Zn,Bd) , np.sqrt(s2y + s2u) )
    return pdf

def pdf_cat(Zn,Bd,s2u,R):
    # Function to compute pdf of an ordinal variable. It returns the whole
    # pdf, a probability vector of length R (number of categories)
    # Input parameters:
    #    Zn: [1*K], feature activation vector
    #    Bd: [K*1], feature weights (dictionary)
    #   s2u: scalar, variance of auxiliary noise
    #     R: number of categories

    pdf = np.zeros(R)
    numMC_samples = 100
    for r in xrange(R):
        tmp = np.ones((1,numMC_samples))
        # we compute the expectation using Monte Carlo samples
        # TODO: check that u does not depend on j
        uV = np.sqrt(s2u) * np.random.randn(numMC_samples) # mean = 0
        for j in xrange(R):
            if (j==r):
                continue
            tmp = tmp * norm.cdf(uV + np.kron( np.ones(numMC_samples), \
                    np.inner(Zn,Bd[:,r]-Bd[:,j]) ) )
                #(Bd(:,r)-Bd(:,j)), 1, numMC_samples) )
        pdf[r] = np.mean(tmp,1)
    return pdf


def pdf_count(X,Zn,Bd,w,s2y, fpos_1_handler):
    pdf = norm.cdf(fpos_1_handler(X+1,w), np.inner(Zn,Bd), np.sqrt(s2y)) \
        - norm.cdf(fpos_1_handler(X,w), np.inner(Zn,Bd), np.sqrt(s2y))
    return pdf

def pdf_ord(Zn,Bd,theta,s2y):
    # Function to compute pdf of an ordinal variable. It returns the whole
    # pdf, a probability vector of length R (number of categories)
    # Input parameters:
    #    Zn: [1*K], feature activation vector
    #    Bd: [K*1], feature weights (dictionary) # TODO: Review dimensions
    # theta: [1*(R-1)]
    #   s2y: scalar, variance of pseudo-observations
    R = len(theta)+1 # number of categories
    pdf = np.zeros(R)
    for r in xrange(R):
        if (r==0):
            a = norm.cdf(theta[0],np.inner(Zn,Bd),np.sqrt(s2y))
            b = 0
        elif (r==(R-1)):
            a = 1
            b = norm.cdf(theta[r-2],np.inner(Zn,Bd),np.sqrt(s2y))
        else:
            a = norm.cdf(theta[r-1],np.inner(Zn,Bd),np.sqrt(s2y))
            b = norm.cdf(theta[r-2],np.inner(Zn,Bd),np.sqrt(s2y))
        pdf[r] = a - b
    return pdf

def pdf_pos(X,Zn,Bd,w,s2y,s2u, fpos_1_handler, dfpos_1_handler):
    # Probability density function for positive real variables
    # Inputs:
    #   X: values in which to evaluate pdf
    #  Zn: [1*K] vector
    #  Bd: [K*1] vector
    #   w: scalar: normalization weight # TODO: Review how to automatize
    # s2y: scalar, variance of pseudo-observations
    # s2u: scalar, variance of auxiliary noise

    pdf = 1.0/np.sqrt(2*np.pi*(s2u+s2y)) * np.exp( \
            -1.0/(2.0*(s2y+s2u)) * (fpos_1_handler(X,w) - np.inner(Zn,Bd))**2)\
            * np.abs( dfpos_1_handler(X,w) )
    return pdf
