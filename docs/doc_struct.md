[**Introduction**](https://ivaleram.github.io/GLFM/) | [**Functions**](doc_functions.html) | [**Data Structures**](doc_struct.html) | [**Demos**](demos.html) | [**FAQ**](FAQ_errors.html)

## Description of data structures defined in GLFM package

Basic definitions
--------------------------
* N: number of observations
* D: number of dimensions
* K: number of inferred latent features
* maxR: maximum number of categories among all categorical variables

Structures (lists in the case of R)
--------------------------
1. **data**: structure with data to analyze/complete
    * X:  N x D observation matrix of N samples and D dimensions
    * C:  1 x D string array indicating type of data for each dimension. There are currently five types of data defined: 'g'= continuous real-valued; 'p'= continuous positive real-valued; 'n'= discrete count data; 'c'= discrete categorical; and 'o'= discrete ordinal
    * ylabel: 1 x D string cell with variable names
    * cat_labels: 1 x D cell with categorical labels (cell is empty for other data types)

2. **hidden**: structure with all latent variables and data-driven parameters
    * Z:  N x K binary matrix of feature assignments (initialization for the IBP)
    * B:  D x K x maxR  weighting tuple
    * theta: D x marR matrix of thresholds for ordinal observations
    * mu: 1 x D vector of shift parameters for internal transformation
    * w: 1 x D vector or scale parameters for internal transformation
    * s2y: 1 x D vector with inferred noise variance for pseudo-observations Y
    * R: 1 x D vector with unique labels for discrete variables (categorical or ordinal data), None otherwise

3. **params**: structure with all simulation parameters (if not defined, the default values listed below will be used):

    * missing = -1; % value for missings
    * alpha = 1; % Concentration parameter for the IBP
    * bias = 0; % number of bias terms, i.e., number of latent features not to be sampled
    * s2u = 0.01; % variance of auxiliary noise
    * s2B = 1; % Variance of the Gaussian prior of the weighing vectors B^d
    * Niter = 1000; % number of Gibbs sampling iterations to run
    * maxK = size(data.X,2); % maximum number of latent features for memory allocation inside C++ routine
    * verbose = 1; % show info in command line

    -- parameters for **optional data-preprocessing** --
    * t = cell(1,size(data.X,2) ); % optional external transformation of obs. X(d) = params.t{d}(Xraw(d)). For example, params.t_1{d} = @(x) log(x + 1);
    * t_1 = cell(1,size(data.X,2) ); % inverse transform. For example, params.t{d} = @(y) exp(y) - 1;
    * dt_1 = cell(1,size(data.X,2) ); % derivative of the inverse transform. For example, params.dt_1{d} = @(x) 1./(x+1);

4. (Only for R) **output**: list containing the output lists **hidden** and **params**.

Examples of data structures
-----------------------------

**Examples of structure "data" in python**

{'cat_labels': [['3', '4'], ['0', '1', '5'], ['alive', 'vascular', 'prostatic ca', 'others'], None, None], 'X': array([[ 3.        ,  0.        ,  1.        ,  2.        ,  0.29998779],
      [ 3.        ,  0.        ,  4.        , 42.        ,  0.69995117],
      [ 3.        ,  5.        ,  2.        ,  3.        ,  0.29998779],
      ...,
      [ 3.        ,  1.        ,  1.        , 11.        ,  0.79992676],
      [ 4.        ,  5.        ,  3.        , 32.        ,  8.3984375 ],
      [ 4.        ,  1.        ,  3.        , 33.        , 22.19921875]]), 'C': 'cocnp', 'ylabel': ['Stage', 'Drug level', 'Prognostic Status', 'Size of Primary Tumor (cm^2)', 'Serum Prostatic Acid Phosphatase']}

{'cat_labels': [['Alaska', 'Alabama', 'Arkansas', 'Arizona', 'California', 'Colorado', 'Connecticut', 'District of Columbia', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Iowa', 'Idaho', 'Illinois', 'Indiana', 'Kansas', 'Kentucky', 'Louisiana', 'Massachusetts', 'Maryland', 'Maine', 'Michigan', 'Minnesota', 'Missouri', 'Mississippi', 'Montana', 'North Carolina', 'North Dakota', 'Nebraska', 'New Hampshire', 'New Jersey', 'New Mexico', 'Nevada', 'New York', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Virginia', 'Vermont', 'Washington', 'Wisconsin', 'West Virginia', 'Wyoming'], None, None, None, None, None, None, None, None, None, None], 'X': array([[2.00000000e+00, 6.10000000e+01, 9.79999971e+00, ...,
        1.23000002e+01, 7.93173981e+01, 2.00017529e+01],
       [2.00000000e+00, 6.70000000e+01, 1.51999998e+01, ...,
        1.65000000e+01, 8.60449753e+01, 1.28612127e+01],
       [2.00000000e+00, 2.90000000e+01, 1.45999999e+01, ...,
        9.80000019e+00, 5.55455017e+01, 4.40413895e+01],
       ...,
       [5.10000000e+01, 9.00000000e+00, 5.29999995e+00, ...,
        2.96000004e+01, 9.77171860e+01, 1.33654103e-01],
       [5.10000000e+01, 4.00000000e+00, 1.38000002e+01, ...,
        2.72000008e+01, 9.37529831e+01, 1.66905105e-01],
       [5.10000000e+01, 3.00000000e+00, 1.27000003e+01, ...,
        2.73999996e+01, 9.80362091e+01, 4.60263900e-02]]), 'C': 'cnppnpppppp', 'ylabel': ['state', 'pop.density', 'age >= 65', 'college', 'income', 'farm', 'democrat', 'republican', 'Perot', 'white', 'black']}

**Examples of structure "params" in python**

{'ext_dataType': [None, None, None, None, 'p'], 'Niter': 100, 'dt_1': [None, None, None, None, <function <lambda> at 0x7fa8e53ce5f0>], 'maxK': 10, 's2B': 1, 'bias': 1, 't': [None, None, None, None, <function <lambda> at 0x7fa8e53ce500>], 's2u': 0.005, 'alpha': 1, 't_1': [None, None, None, None, <function <lambda> at 0x7fa8e53ce578>]}

{'s2B': 1, 'Niter': 100, 'dt_1': [None, <function <lambda> at 0x7f345f9015f0>, None, None, <function <lambda> at 0x7f345f901758>, <function <lambda> at 0x7f345f9018c0>, None, None, None, <function <lambda> at 0x7f345f901b90>, <function <lambda> at 0x7f345f901a28>], 'maxK': 10, 'bias': 1, 't': [None, <function <lambda> at 0x7f345f901500>, None, None, <function <lambda> at 0x7f345f9016e0>, <function <lambda> at 0x7f345f901848>, None, None, None, <function <lambda> at 0x7f345f901b18>, <function <lambda> at 0x7f345f9019b0>], 's2u': 0.005, 'alpha': 1, 't_1': [None, <function <lambda> at 0x7f345f901578>, None, None, <function <lambda> at 0x7f345f901668>, <function <lambda> at 0x7f345f9017d0>, None, None, None, <function <lambda> at 0x7f345f901aa0>, <function <lambda> at 0x7f345f901938>], 'ext_dataType': [None, 'p', None, None, 'p', 'p', None, None, None, 'p', 'p']}

**Examples of structure "hidden" in python**

{'Z': array([[1., 0.],
       [1., 0.],
       [1., 0.],
       ...,
       [1., 1.],
       [1., 0.],
       [1., 0.]])}


