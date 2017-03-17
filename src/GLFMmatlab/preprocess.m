function [Xnorm,suffStats] = preprocess(X,C)
    % Function to compute normalization mapping functions
    % Inputs:
    %   X: Raw matrix [N*D]
    %   C: String [1*D] with datatypes
    %
    % Outputs:
    %   Xnorm: Normalized matrix
    %   suffStats: Cell [1*D] with array of parameters for each datatype
    %              transformation
    %
    W = zeros(1,D); % vector of weights for transformation
    for d=1:D
        % normalize data (such that it occupies interval (0,max(X(d))]
        if ((data.C(d) == 'n') || (data.C(d) == 'c') || (data.C(d) == 'o'))   && (min(data.X(:,d)) > 1)
            offset = min(Xmiss(Xmiss(:,d) ~= missing,d));
            Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 1;
        end
        if (data.C(d) == 'p') && (min(data.X(:,d)) > 0)
            offset = min(data.X(:,d));
            Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 10^-6;
        end
        W(d) = 2/max( Xmiss(Xmiss(:,d) ~= missing,d) );
        
        if ((data.C(d) == 'n' || (data.C(d) == 'c') || (data.C(d) == 'o')) && (min(data.X(:,d)) == 0))
            Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 1;
        elseif (data.C(d) == 'p') && (min(data.X(:,d)) == 0)
            Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 10^-6;
        end
    end