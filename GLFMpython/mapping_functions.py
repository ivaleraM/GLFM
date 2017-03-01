# -------------------------------------------------------------------------
# Module with all available mapping functions necessary for the Generalized
# Latent Feature Model (GLFM)
# -------------------------------------------------------------------------

# --------------------
# Auxiliary Functions
# --------------------

def exp_count(Zn,Bd, params):
    # function to compute expectation of random variable

    s2y = params{1}
    fpos_1_handler = params{2}
    w = params{3} # function to compute normalization weights w
    maxX = params{4}

    xin = 1:1:maxX

    func = @(xin) pdf_count(xin,Zn,Bd,w, s2y, fpos_1_handler)

    expect = xin * func(xin) # transpose '
    return expect

def rnd_pos(Zn,Bd,numSamples,params):
    # function to get random samples from the distribution

    s2y = params{1}
    s2u = params{2}
    fpos_1_handler = params{3}
    dfpos_1_handler = params{4}
    w = params{5} # function to compute normalization weights w
    maxX = params{6}

    func = @(xin) log( pdf_pos(xin,Zn,Bd,w, s2y ,s2u, fpos_1_handler, dfpos_1_handler) )
    a = 10^-6
    b = maxX
    domain = [a,b+5]

    samples = ars(func, a, b, domain, numSamples, [])
    samples = exp(samples)
    return samples

def rnd_real(Zn,Bd,numSamples, s2y, s2u):
    # function to get random samples from the distribution
    x = sqrt(s2u+s2y) .* randn(1,numSamples) + Zn * Bd
    return x

# --------------------
# Mapping Functions
# --------------------

def freal(y, s2u):
    # Mapping function for real values
    x = y + sqrt(s2u) * randn(size(y))
    return x

def fcat(y):
    # input argument y: [N*R]
    # output: x [N*1]
    [mm, x] = max(y,[],2)
    return x

def fcount(y,w):
    # Inputs:
    #   x: [N*D]
    x = floor( fpos(y, w) )
    return x

def ford(y, theta):
    # Mapping function for ordinal data
    # Inputs:
    #       y: [1*R] Pseudo-observations
    #   theta: [1*(R-1)] Thresholds that divide the real line into R regions
    for r in xrange(len(theta)): #= 1:length(theta)
        val = theta()
        if (y < val):
            x = r
            break
        if r == size(theta,3):
            x = r+1
    return x

def fpos(y, w):
    # Inputs:
    #   y: 
    #   w: scalar

    ##x = log( exp( w * y ) + 1)
    #x =  y^2 /w
    x = log(exp(y)+1) / w
    return x

def fpos_1(x,w):
    # w: scalar

    #y = sqrt(w*x)
    # y = log( exp( x ) + 1) / w
    y = log( exp( w*x ) -1 )
    return y

def fpos_1_xi(x,w):
    # w: scalar

    y = sqrt(w*x)
    # y = log( exp( x ) + 1) / w
    # y = log( exp( w*x ) -1 )
    return y

def fpos_xi(y, w):
    # Inputs:
    #   y: 
    #   w: scalar
    # TODO: make weights W input-dependent
    #x = log( exp( w * y ) + 1)

    x =  y^2 /w
    #x = log(exp(y)+1) / w
    return x

def dfpos_1(x, w):
    # w: scalar
    #y = 1./sqrt(x)
    #y = -0.5*w^0.5 * x.^(-0.5)
    y = w ./ ( 1 - exp( w .* x) )
    return y

def dfpos_1_xi(x, w):
    # w: scalar

    #y = 1./sqrt(x)
    y = -0.5*w^0.5 * x.^(-0.5)
    #y = w ./ ( 1 - exp( w .* x) )
    return y

# ------------------------------------------
# Functions to compute pdf values
# ------------------------------------------

def pdf_real(X, Zn,Bd,s2y,s2u):
    # Probability Density Function for real variables
    pdf = normpdf(X, Zn* Bd, sqrt(s2y + s2u) )
    return pdf

def pdf_cat(Zn,Bd,s2u,R):
    # Function to compute pdf of an ordinal variable. It returns the whole
    # pdf, a probability vector of length R (number of categories)
    # Input parameters:
    #    Zn: [1*K], feature activation vector
    #    Bd: [K*1], feature weights (dictionary)
    #   s2u: scalar, variance of auxiliary noise
    #     R: number of categories

    pdf = zeros(1,R)
    numMC_samples = 100
    # TODO: perform pdf computation in log space ?
    for r=1:R
        tmp = ones(1,numMC_samples)
        # we compute the expectation using Monte Carlo samples
        uV = sqrt(s2u) * randn(1,numMC_samples) # mean = 0 # TODO: check that u does not depend on j
        for j=1:R
            if (j==r)
                continue
            end
            tmp = tmp .* normcdf(uV + repmat( Zn * (Bd(:,r)-Bd(:,j)), 1, numMC_samples) )
        end
        pdf(r) = mean(tmp,2)
    end
    return pdf


def pdf_count(X,Zn,Bd,w,s2y, fpos_1_handler):

pdf = normcdf(fpos_1_handler(X+1,w), Zn*Bd, sqrt(s2y)) - ...
    normcdf(fpos_1_handler(X,w), Zn*Bd, sqrt(s2y))
    return pdf

def pdf_ord(Zn,Bd,theta,s2y):
    # Function to compute pdf of an ordinal variable. It returns the whole
    # pdf, a probability vector of length R (number of categories)
    # Input parameters:
    #    Zn: [1*K], feature activation vector
    #    Bd: [K*1], feature weights (dictionary)
    # theta: [1*(R-1)]
    #   s2y: scalar, variance of pseudo-observations
    R = length(theta)+1 # number of categories
    pdf = zeros(1,R)
    for r=1:R
        if (r==1)
            a = normcdf(theta(1),Zn * Bd,s2y)
            b = 0
        elseif (r== R)
            a = 1
            b = normcdf(theta(r-1),Zn * Bd,s2y)
        else
            a = normcdf(theta(r),Zn * Bd,s2y)
            b = normcdf(theta(r-1),Zn * Bd,s2y)
        end
        pdf(r) = a - b
    end
end
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

    pdf = 1/sqrt(2*pi*(s2u+s2y)) .* exp( -1/(2*(s2y+s2u)) * (fpos_1_handler(X,w) - Zn*Bd).^2 ) ...
        .* abs( dfpos_1_handler(X,w) )
    return pdf
