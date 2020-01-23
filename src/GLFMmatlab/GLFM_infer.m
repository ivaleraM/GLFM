function [hidden,data_transformed] = GLFM_infer(data,varargin)
    % Wrapper .m function to call .cpp MATLAB wrapper (simplifies call)
    % Three possible calls:
    %           hidden = GLFM_infer(data)
    %           hidden = GLFM_infer(data,hidden)
    %           hidden = GLFM_infer(data,hidden,params)
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
    %       params: structure with simulation parameters and hyperparameters
    %   Outputs:
    %       hidden: hidden structure with samples from posterior, as well as
    %               inner transformation parameters
    %           hidden.Z    : N*K feature activation matrix
    %           hidden.B    : K*D feature effect matrix
    %           hidden.theta: D*maxR matrix with aux. variables for ordinal data
    %           hidden.mu   : 1*D shift parameter for inner transformation
    %           hidden.w    : 1*D scale parameter for inner transformation
    %           hidden.s2y  : 1*D per-dim inferred noise variance for pseudo-obs

    switch length(varargin)
        case 0
            hidden = [];
            % initialize default values for params
            params = init_default_params(data, []);

        case 1, hidden = varargin{1};
            params = init_default_params(data, []);

        case 2, hidden = varargin{1}; params = varargin{2};
            params = init_default_params(data, params); % eventually complete params structure

        otherwise
            error('Incorrect number of input parameters: should be 0, 1 or 2');
    end

    D = size(data.X,2);
    N = size(data.X,1);

    % initialize Z if necessary
    if isempty(hidden)
        clearvars hidden
        hidden.Z = double(rand(N,2)>0.8);
        if (params.bias == 1)
            hidden.Z = [ones(N,1) hidden.Z];
        elseif (params.bias > 1)
            error('There is more than 1 bias specified, but structure Z has not been initialized...');
        end
    end

    % replace missings
    data_transformed = data;
    data_transformed.X(isnan(data_transformed.X)) = params.missing;

    % change labels for categorical and ordinal vars such that > 0
    %V_offset = zeros(1,D);
    for d=1:D
        if (data_transformed.C(d) == 'c') || (data_transformed.C(d) == 'o')
            mask = data_transformed.X(:,d) ~= params.missing;
            %V_offset(d) = min( data_transformed.X(mask,d) );
            %data_transformed.X(mask,d) = data_transformed.X(mask,d) - V_offset(d) + 1;
            uniqueVal= unique(data_transformed.X(mask,d));
            Xaux=[];
            for i=1:length(uniqueVal)
                Xaux(data_transformed.X(:,d)==uniqueVal(i)) = i;
            end
            Xaux(~mask)=params.missing;
            data_transformed.X(:,d)= Xaux;
        end
    end

    % eventually, apply external transform
    for r=1:size(data_transformed.X,2)
        if ~isempty(params.t{r})
            data_transformed.X(:,r) = params.t_1{r}(data_transformed.X(:,r)); % work in logarithm space better
            data_transformed.C(r) = params.ext_dataType{r};
        end
    end

    func = 1*ones(1,size(data_transformed.X,2)); % type of internal transformation for positive real-valued data

    %% call .cpp wrapper function
    tic;
    [Z B theta mu w s2Y]= IBPsampler(data_transformed.X,data_transformed.C, hidden.Z, params.bias, func, ...
        params.s2u, params.s2B, params.alpha, params.Niter, params.maxK, params.missing);
    hidden.time = toc;
    fprintf('Elapsed time %.2f seconds.', hidden.time );

    %% prepare output structure
    hidden.Z = Z'; % it returns a K*N matrix, should be inverted
    hidden.B = B;
    hidden.theta = theta;
    hidden.mu = mu; % mapping functions info (necessary for plotting)
    hidden.w = w;
    hidden.s2Y = s2Y;
    hidden.R = ones(1,D);
    for d=1:D
        if (data_transformed.C(d) == 'c') || (data_transformed.C(d) == 'o')
            hidden.R(d) = length(unique(data_transformed.X(data_transformed.X(:,d)~= params.missing,d) ));
        end
    end
    if params.verbose
        hidden.time
        sum(hidden.Z)
    end
end
