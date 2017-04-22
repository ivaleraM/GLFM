function hidden = IBPsampler_run(data,varargin)
    % Wrapper .m function to call .cpp MATLAB wrapper (simplifies call)
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
    
    switch length(varargin)
        case 0
            hidden = [];
            % default values for params
            params = init_default_params(data, []);
            
        case 1, hidden = varargin{1};
            params = init_default_params(data, []);
            
        case 2, hidden = varargin{1}; params = varargin{2};
            params = init_default_params(data, params); % eventually complete params structure
            
        otherwise
            error('Incorrect number of input parameters: should be 1, 2 or 3');
    end
    
    D = size(data.X,2);
    N = size(data.X,1);
    
    % initialize Z
    if isempty(hidden)
        clearvars hidden
        hidden.Z = double(rand(N,2)>0.8);
        if (params.bias == 1)
            hidden.Z = [ones(N,1) hidden.Z];
        elseif (params.bias > 1)
            error('There is more than 1 bias specified, but structure Z has not been initialized...');
        end
    end
    
    % replace missings + preprocess
    data.X(isnan(data.X)) = params.missing;
    
    % change labels for categorical and ordinal vars such that > 0
    V_offset = zeros(1,D);
    for d=1:D
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            mask = data.X(:,d) ~= params.missing;
            V_offset(d) = min( data.X(mask,d) );
            data.X(mask,d) = data.X(mask,d) - V_offset(d) + 1;
        end
    end
    
    % eventually, apply external transform
    for r=1:size(data.X,2)
        if ~isempty(params.t{r})
            data.X(:,r) = params.t_1{r}(data.X(:,r)); % work in logarithm space better
            data.C(r) = params.ext_dataType{r};
        end
    end
    
    %% call .cpp wrapper function
    tic;
    [Z B theta mu w s2Y]= IBPsampler(data.X,data.C, hidden.Z, params.bias,params.func, params.s2Y, ...
        params.s2u, params.s2B, params.alpha, params.Niter, params.maxK, params.missing);
    hidden.time = toc;
    
    hidden.Z = Z'; % it returns a K*N matrix, should be inverted
    hidden.B = B;
    hidden.theta = theta;
    hidden.mu = mu; % mapping functions info (necessary for plotting)
    hidden.w = w;
    hidden.s2Y = s2Y;
    hidden.R = ones(1,D);
    for d=1:D
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            hidden.R(d) = length(unique(data.X(data.X(:,d)~= params.missing,d) ));
        end
    end
    if params.verbose
        hidden.time
        sum(hidden.Z)
    end
end