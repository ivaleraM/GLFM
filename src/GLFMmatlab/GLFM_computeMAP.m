function X_map = GLFM_computeMAP(C, Zp, hidden, params, varargin)
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
    % ----------------(optional) ------------------
    %       - idxsD: dimensions to infer
    %
    % Outputs:
    %   X_map: P*Di matrix with MAP estimate where Di = length(idxsD)

    if (length(varargin) == 1)
        idxsD = varargin{1};
    elseif (length(varargin) > 1)
        error('Too many input arguments');
    else
        idxsD = 1:size(hidden.B,1);
    end
    P = size(Zp,1);
    K = size(hidden.B,2);
    if (size(Zp,2) ~= K)
        error('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different');
    end

    X_map = zeros(P,length(idxsD)); % output
    for dd=1:length(idxsD) % for each dimension
        d = idxsD(dd);
        if isfield(params,'t') % if external transformations have been defined
            if ~isempty(params.t{d}) % there is an external transform for data type d
                C(d) = params.ext_dataType{d}; % set new type of data
            end
        end
        switch C(d)
            case 'g', X_map(:,dd) = f_g( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
            case 'p', X_map(:,dd) = f_p( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
            case 'n', X_map(:,dd) = f_n( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
            case 'c', X_map(:,dd) = f_c( Zp * squeeze(hidden.B(d,:,1:hidden.R(d))) );
            case 'o', X_map(:,dd) = f_o( Zp * squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)) );
            otherwise
                error('Unknown data type');
        end
        if (sum(isnan(X_map(:,dd))) > 0)
            warning('Some values are nan!');
        end
        if isfield(params,'t')
            if ~isempty(params.t{d})
                X_map(:,dd) = params.t{d}( X_map(:,dd) );
            end
        end
    end
