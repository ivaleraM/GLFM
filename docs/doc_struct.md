[MATLAB functions](doc_matlab.html) | [PYTHON functions](doc_python.html) | [Demos](demos.html) | [FAQ](FAQ_errors.html)

Description of all used structures. See full description of [available functions](doc_matlab.html) in MATLAB.

N: number of observations
D: number of dimensions
K: number of inferred latent features
maxR: maximum number of categories among all categorical variables

1. **data**: structure with data to analyze/complete
    * X: ´N x D´ observation matrix of N samples and D dimensions
    * C: ´1 x D´ string array indicating type of data for each dimension
    * ylabel: ´1 x D´ string cell with variable names
    * cat_labels: ´1 x D´ cell with categorical labels (cell is empty for other data types)

2. **hidden**: structure with all latent variables and data-driven parameters
    * Z: ´N x K´ binary matrix of feature assignments (initialization for the IBP)
    * B: ´D x K x maxR´ feature effect matrix
    * theta: ´D x marR´ matrix of auxiliary variables
    * mu: ´1 x D´ vector of shift parameters for internal transformation
    * w: ´1 x D´ vector or scale parameters for internal transformation
    * s2y: ´1 x D´ vector with inferred noise variance for pseudo-observations Y

3. **params**: structure with all simulation parameters (if not defined, the default values listed below will be used):

        missing = -1; end % value for missings
        alpha = 1; end % Concentration parameter for the IBP
        bias = 0; end % number of latent features that should not be sampled
        s2Y = 0.5; end % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
        s2u = 0.01; end % variance of auxiliary noise
        s2B = 1; end % Variance of the Gaussian prior of the weigting matrices B
        Niter = 1000; end % number of iterations to run
        maxK = size(data.X,2); end % max number of latent features for memory allocation inside C++ routine
        verbose = 1; end % plot info in command line

        eters for optional external transformation
        t = cell(1,size(data.X,2) ); end % eventual external transform of obs. X = params.t{d}(Xraw)
        t_1 = cell(1,size(data.X,2) ); end % inverse transform
        dt_1 = cell(1,size(data.X,2) ); end % derivative of the inverse transform
