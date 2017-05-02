
## Documentation for PYTHON functions

**def infer**(Xin,Cin,Zin,bias=0, s2u=0.001, s2B=1.0,
        alpha=1.0, Nsim=100, maxK=50, missing=-1, verbose=0):

    """
    Python wrapper to launch inference routine for the GLFM model.
    Inputs:
        Xin: input data D*N where:
                N = number of observations
                D = number of dimensions
        Cin: array indicating types of data ('g': real,'p': positive real-valued,
            'c': categorical; 'o': ordinal; 'n': count data)
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
