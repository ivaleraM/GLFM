[**Introduction**](https://ivaleram.github.io/GLFM/) | [**Functions**](doc_functions.html) | [**Data Structures**](doc_struct.html) | [**Demos**](demos.html) | [**FAQ**](FAQ_errors.html)

## Documentation for MATLAB functions

See complete description of [data structures](doc_struct.html).

**function hidden = GLFM.infer(data,varargin)**

    % Wrapper .m function to call .cpp MATLAB wrapper (simplifies call)
    % Three possible calls:
    %           hidden = GLFM.infer(data)
    %           hidden = GLFM.infer(data,hidden)
    %           hidden = GLFM.infer(data,hidden,params)
    %
    %   Inputs:
    %       data: structure with all input data information
    %           (*) data.X: N*D observation matrix (raw)
    %           (*) data.C: 1*D string array with input data types
    %       (*) mandatory
    %
    %       ------------- optional ---------------------
    %
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with sim. parameters and hyperparameters

**function [Xcomplete,hidden] = GLFM.complete(data,varargin)**

    % Function to complete a matrix that has missing values
    % Possible calls:
    %           [Xcomplete,hidden]  = GLFM.complete(data)
    %           [Xcomplete,hidden]  = GLFM.complete(data,hidden) % init hidden.Z externaly
    %           [Xcomplete,hidden]  = GLFM.complete(data,[],params) % struc. with parameters
    %           [Xcomplete,hidden]  = GLFM.complete(data,hidden,params)
    %
    %   Inputs:
    %       data: structure with all input data information
    %           (*) data.X: NxD observation matrix (raw) with missings
    %           (*) data.C: 1xD string array with input data types
    %       (*) mandatory
    %
    %       ------------- optional ---------------------
    %
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with sim. parameters and hyperparameters
    %   Output:
    %       Xmap: NxD input matrix with imputed missing values
    %       hidden: structure with latent parameters (same output as
    %       IBPsampler_infer function).

**function X_map = GLFM.computeMAP(C, Zp, hidden, params)**

    % Function to generate the MAP solution corresponding to patterns in Zp
    % Inputs:
    %   C: 1*D string with data types, D = number of dimensions
    %   Zp: P * K matrix of feature activation for which to compute the MAP estimate
    %       (P is the number of obs.)
    %   hidden: structure with latent variables learned by the model
    %       - B: latent feature matrix (D * K * maxR)  where
    %               D: number of dimensions
    %               K: number of latent variables
    %            maxR: maximum number of categories across all dimensions
    %       - mu: 1*D shift parameter
    %       - w:  1*D scale parameter
    %       - theta: D*maxR matrix of auxiliary vars (for ordinal variables)
    %
    % Outputs:
    %   X_map: P*D matrix with MAP estimate


**function [xd, pdf] = IBPsampler_computePDF(data, Zp, hidden, params, d)**

    % Function to generate the PDF solutions corresponding to patterns in
    % Zp, and dimension d
    % Inputs:
    %   data.X: N*D data matrix (necessary to compute domain of x
    %   data.C: 1*D string with data types, D = number of dimensions
    %   Zp: P * K matrix of patterns (P is the number of patterns)
    %   hidden.B: latent feature matrix (D * K * maxR)   
    %   hidden.mu: 1*D shift parameter for internal transformation
    %   hidden.w:  1*D scale parameter for internal transformation
    %   hidden.Theta: D*maxR matrix of aux. variables for ordinal data
    %                 where R is the max number of categories
    %  in params, we can specify number of points numS to compute (for 'g' and 'p' types)
    %
    % Outputs:
    %   xd: 1*numS where numS is the number of points to compute
    %  pdf: P*numS where P is the number of patterns to consider
