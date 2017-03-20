function hidden = run_IBPsampler(data,varargin)
    % Wrapper .m function to call .cpp functions and make call simpler
    %
    %   Inputs:
    %       data: structure with all input data information
    %           (*) data.X: N*D observation matrix (raw)
    %           (*) data.C: 1*D string array with input data types
    %
    %       ------------- optional ---------------------
    %
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with sim. parameters and hyperparameters
    
    switch length(varargin)
        case 0
            % initialize Z
            
            % default values for params
            params.bias = 0;
            params.s2Y = 0.1;
            params.s2u = 0.1;
            params.s2B = 1;
            params.alpha = 1;
            params.Niter = 100;
            params.maxK = size(data.X,2);
            params.missing = -1;
        case 1, hidden = varargin{1};
        case 2, hidden = varargin{1}; params = varargin{2};
        otherwise
            error('Incorrect number of input parameters: should be 1, 2 or 3');
    end
    
    % replace missings + preprocess
    data.X(isnan(data.X)) = params.missing;
    
    %[Xnorm,suffStats] = preprocess(X,C)
    [Xmiss, suffStats] = preprocess(data.X, data.C, params.missing);
    
    % call .cpp wrapper function
    [Z B Theta]= IBPsampler(Xmiss,data.C, hidden.Z, params.bias, params.s2Y, ...
        params.s2u, params.s2B, params.alpha, params.Niter, params.maxK, params. missing);
    
    hidden.Z = Z'; % it returns a K*N matrix, should be inverted
    hidden.B = B;
    hidden.Theta = Theta;
    hidden.suffStats = suffStats; % mapping functions info (necessary for plotting)
    
end