[**Introduction**](https://ivaleram.github.io/GLFM/) | [**Functions**](doc_functions.html) | [**Data Structures**](doc_struct.html) | [**Demos**](demos.html) | [**FAQ**](FAQ_errors.html)

## Description of the data structures defined in GLFM package 

Basic definitions
--------------------------
* N: number of observations
* D: number of dimensions
* K: number of inferred latent features
* maxR: maximum number of categories among all categorical variables

Structures
--------------------------
1. **data**: structure with data to analyze/complete
    * X:  N x D observation matrix of N samples and D dimensions
    * C:  1 x D string array indicating type of data for each dimension
    * ylabel: 1 x D string cell with variable names
    * cat_labels: 1 x D cell with categorical labels (cell is empty for other data types)

2. **hidden**: structure with all latent variables and data-driven parameters
    * Z:  N x K binary matrix of feature assignments (initialization for the IBP)
    * B:  D x K x maxR  weighting tuple
    * theta: D x marR matrix of thresholds for ordinal observations
    * mu: 1 x D vector of shift parameters for internal transformation
    * w: 1 x D vector or scale parameters for internal transformation
    * s2y: 1 x D vector with inferred noise variance for pseudo-observations Y

3. **params**: structure with all simulation parameters (if not defined, the default values listed below will be used):

        missing = -1; % value for missings
        alpha = 1; % Concentration parameter for the IBP
        bias = 0; % number of bias terms, i.e., number of latent features not to be sampled
        s2Y = 0.5; % Variance of the Gaussian prior on the auxiliary variables (pseudo-observations) Y
        s2u = 0.01; % variance of auxiliary noise
        s2B = 1; % Variance of the Gaussian prior of the weighing vectors B^d
        Niter = 1000; % number of Gibbs sampling iterations to run
        maxK = size(data.X,2); % maximum number of latent features for memory allocation inside C++ routine
        verbose = 1; % show info in command line

        * parameters for optional data-preprocessing
        t = cell(1,size(data.X,2) ); end % eventual external transform of obs. X = params.t{d}(Xraw)
        t_1 = cell(1,size(data.X,2) ); end % inverse transform
        dt_1 = cell(1,size(data.X,2) ); end % derivative of the inverse transform
