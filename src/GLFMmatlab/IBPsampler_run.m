function hidden = IBPsampler_run(data,varargin)
    % Wrapper .m function to call .cpp MATLAB wrapper (simplifies call)
    %
    %   Inputs:
    %       data: structure with all input data information
    %           (*) data.X: N*D observation matrix (raw)
    %           (*) data.C: 1*D string array with input data types
    %       (*) means mandatory
    %
    %       ------------- optional ---------------------
    %
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with sim. parameters and hyperparameters
    
    switch length(varargin)
        case 0
            hidden.Z = [];
            % default values for params
            params = init_default_params(data, []);
            
        case 1, hidden = varargin{1};
            params = init_default_params(data, []);
            
        case 2, hidden = varargin{1}; params = varargin{2};
            params = init_default_params(data, params); % eventually complete params structure
            
        otherwise
            error('Incorrect number of input parameters: should be 1, 2 or 3');
    end
    
    % initialize Z
    if isempty(hidden.Z)
        hidden.Z = double(rand(N,2)>0.8);
        if (params.bias == 1)
            hidden.Z = [ones(N,1) hidden.Z];
        elseif (params.bias > 1)
            error('There is more than 1 bias specified, but structure Z has not been initialized...');
        end
    end
    
    % replace missings + preprocess
    data.X(isnan(data.X)) = params.missing;
    
    %[Xmiss, suffStats] = preprocess(data.X, data.C, params.missing);
    
    %% call .cpp wrapper function
    tic;
    [Z B theta mu w]= IBPsampler(data.X,data.C, hidden.Z, params.bias, ones(1,size(data.X,2)), params.s2Y, ...
        params.s2u, params.s2B, params.alpha, params.Niter, params.maxK, params.missing);
    hidden.time = toc;
    
    hidden.Z = Z'; % it returns a K*N matrix, should be inverted
    hidden.B = B;
    hidden.theta = theta;
    hidden.mu = mu; % mapping functions info (necessary for plotting)
    hidden.w = w;
    
    if params.verbose
        hidden.time
        sum(hidden.Z)
    end
end