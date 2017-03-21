function params = init_default_params(data, params)
    % Function to initialize or complete params structure with all
    % simulation parameters and hyperparameters for the GLFM
    % data is an input parameter in order to be able to initialize
    % params with data-based values
    
    if ~isfield(params,'missing'), params.missing = -1; end
    if ~isfield(params,'alpha'), params.alpha = 1; end % Concentration parameter for the IBP
    if ~isfield(params,'bias'), params.bias = 0; end % number of latent features that should not be sampled
    if ~isfield(params,'s2Y'), params.s2Y = 0.5; end % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
    if ~isfield(params,'s2u'), params.s2u = 0.01; end % variance of auxiliary noise
    if ~isfield(params,'s2B'), params.s2B = 1; end % Variance of the Gaussian prior of the weigting matrices B
    if ~isfield(params,'Niter'), params.Niter = 1000; end % number of iterations to run
    if ~isfield(params,'maxK'), params.maxK = size(data.X,2); end % max number of latent features for memory allocation inside C++ routine
