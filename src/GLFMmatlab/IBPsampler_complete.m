function [Xcompl,hidden] = IBPsampler_complete(data,varargin)
    % Function to complete a matrix that has missing values
    % Possible calls:
    %           hidden = IBPsampler_infer(data)
    %           hidden = IBPsampler_infer(data,hidden) % init hidden.Z externaly
    %           hidden = IBPsampler_infer(data,[],params) % struc. with parameters
    %           hidden = IBPsampler_infer(data,hidden,params)
    %
    %   Inputs:
    %       data: structure with all input data information
    %           (*) data.X: N*D observation matrix (raw) with missings
    %           (*) data.C: 1*D string array with input data types
    %       (*) mandatory
    %
    %       ------------- optional ---------------------
    %
    %       hidden: hidden structure to initialize inference algorithm
    %             hidden.Z: feature assignment N*K matrix
    %       params: structure with sim. parameters and hyperparameters
    %   Output:
    %       Xcompl: N*D input matrix with imputed missing values
    %       hidden: structure with latent parameters (same output as
    %       IBPsampler_infer function).
    
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
            error('Incorrect number of input parameters: should be 1, 2 or 3');
    end
    
    % Run inference
    hidden = IBPsampler_infer(data,hidden,params);
    
    % just in case there is any missing
    data.X(isnan(data.X)) = params.missing;
    
    [xx_miss, yy_miss] = find(data.X == params.missing);
    
    Xcompl = data.X;
    for i=1:length(xx_miss) % for each missing value, compute MAP estimate
        Xcompl(xx_miss(i),yy_miss(i)) = IBPsampler_computeMAP( data.C, hidden.Z(xx_miss(i),:), hidden, params, yy_miss(i));
    end
end
