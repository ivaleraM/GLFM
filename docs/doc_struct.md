Description of all used structures in MATLAB

N: number of observations
D: number of dimensions
K: number of inferred latent features
maxR: maximum number of categories among all categorical variables

1. **data**: structure with data to analyze/complete
    * data.X: ´N x D´ observation matrix of N samples and D dimensions
    * data.C: ´1 x D´ string array indicating type of data for each dimension
    * data.ylabel: ´1 x D´ string cell with variable names
    * data.cat_labels: ´1 x D´ cell with categorical labels (cell is empty for other data types)

2. **hidden**: structure with all latent variables and data-driven parameters
    * hidden.Z: ´N x K´ binary matrix of feature assignments (initialization for the IBP)
    * hidden.B: ´D x K x maxR´ feature effect matrix
    * hidden.theta: ´D x marR´ matrix of auxiliary variables
    * hidden.mu: ´1 x D´ vector of shift parameters for internal transformation
    * hidden.w: ´1 x D´ vector or scale parameters for internal transformation
    * hidden.s2y: ´1 x D´ vector with inferred noise variance for pseudo-observations Y

3. **params**: structure with all simulation parameters (if not defined, the default values listed below will be used):

        params.missing = -1; end % value for missings
        params.alpha = 1; end % Concentration parameter for the IBP
        params.bias = 0; end % number of latent features that should not be sampled
        params.s2Y = 0.5; end % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
        params.s2u = 0.01; end % variance of auxiliary noise
        params.s2B = 1; end % Variance of the Gaussian prior of the weigting matrices B
        params.Niter = 1000; end % number of iterations to run
        params.maxK = size(data.X,2); end % max number of latent features for memory allocation inside C++ routine
        params.verbose = 1; end % plot info in command line

        % parameters for optional external transformation
        params.t = cell(1,size(data.X,2) ); end % eventual external transform of obs. X = params.t{d}(Xraw)
        params.t_1 = cell(1,size(data.X,2) ); end % inverse transform
        params.dt_1 = cell(1,size(data.X,2) ); end % derivative of the inverse transform
