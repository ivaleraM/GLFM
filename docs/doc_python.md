[MATLAB functions](doc_matlab.html) | [PYTHON functions](doc_python.html) | [Demos](demos.html) | [FAQ](FAQ_errors.html)

## Documentation for PYTHON functions

**def infer(Xin,Cin,Zin,bias=0, s2u=0.001, s2B=1.0,
        alpha=1.0, Nsim=100, maxK=50, missing=-1, verbose=0):**

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

**def complete(Xmiss, C, bias=0, s2Y=1, s2u=1, s2B=1, alpha=1, Niter=50, missing=-1):**

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

**def computeMAP(C, Zp, hidden, params, idxsD=[]):**

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

