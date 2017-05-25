[**Introduction**](https://ivaleram.github.io/GLFM/) | [**Functions**](doc_functions.html) | [**Data Structures**](doc_struct.html) | [**Demos**](demos.html) | [**FAQ**](FAQ_errors.html)

## Description of the functions available in the GLFM package

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
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with simulation parameters and hyperparameters
    % Output:
    %    hidden:  structure with latent variables learned by the model
    %	          (same output as GLFM.infer function).


**function [Xcomplete,hidden] = GLFM.complete(data,varargin)**

    % Function to complete a matrix that has missing values with their MAP solution
    % Possible calls:
    %           [Xcomplete,hidden]  = GLFM.complete(data)
    %           [Xcomplete,hidden]  = GLFM.complete(data,hidden) % init hidden.Z externaly
    %           [Xcomplete,hidden]  = GLFM.complete(data,[],params) % struc. with parameters
    %           [Xcomplete,hidden]  = GLFM.complete(data,hidden,params)
    %
    % Inputs:
    %    data: structure with all input data information
    %          (*) data.X: NxD observation matrix (raw) with missings
    %          (*) data.C: 1xD string array with input data types
    %    (*) mandatory
    %
    %       ------------- optional ---------------------
    %    hidden: hidden structure to initialize inference algorithm
    %            hidden.Z: feature assignment N*K matrix
    %    params: structure with simulation parameters and hyperparameters
    % Output:
    %    Xcomplete: NxD input matrix with imputed missing values
    %    hidden:  structure with latent variables learned by the model
    %	          (same output as GLFM.infer function).

**function X_map = GLFM.computeMAP(C, Zp, hidden, params)**

    % Function to generate the MAP solution corresponding to patterns in Zp
    % Inputs:
    %   C: 1*D string with data types, D = number of dimensions
    %   Zp: P * K matrix of feature activation for which to compute the MAP estimate
    %       (P is the number of obs.)
    %   hidden: structure with latent variables learned by the model
    %   params: structure of simulation parameters and hyperparameters
    %
    % Outputs:
    %   X_map: P*D matrix with MAP estimate


**function [xd, pdf] = GLFM.computePDF(data, Zp, hidden, params, d)**

    % Function to generate the PDF solutions corresponding to patterns in
    % Zp, and dimension d
    % Inputs:
    %   data: data structure
    %   Zp: P * K matrix of patterns (P is the number of patterns)
    %   hidden: structure with latent variables learned by the model
    %   params: structure of simulation parameters and hyperparameters
    %   patterns: numP*K list of patterns to plot
    %
    % Outputs:
    %   xd: 1*numS where numS is the number of points to compute
    %   pdf: P*numS where P is the number of patterns to consider
 
**function [] = GLFM.plotPatterns(data, hidden, params, Zp, varargin)**

    % Function to plot the inferred distribution and empirical histogram for each dimension of the observations.
    % Inputs:
    %   data: data structure
    %   hidden: structure with latent variables learned by the model
    %   params: structure of simulation parameters and hyperparameters
    %   Zp: P * K matrix of patterns (P is the number of patterns)
    %   ------ (optional) ------
    %   colors: list of colors to plot
    %   styles: list of styles for each line (for plot, not bar)
    %   leg: legend to use (by default, use patterns as legend)
    %   idxD: array of dimensions to plot
    %
    % Outputs:
    %   void
