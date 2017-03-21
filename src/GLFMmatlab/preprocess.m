function [X,suffStats] = preprocess(X,C,missing)
    % Function to compute normalization previous to mapping functions
    % Inputs:
    %   X: Raw matrix [N*D]
    %   C: String [1*D] with datatypes
    %
    % Outputs:
    %   Xnorm: Normalized matrix
    %   suffStats: Cell [1*D] with array of parameters for each datatype
    %              transformation
    
    D = size(X,2);
    suffStats = cell(1,D);
    for d=1:D
        mask = X(:,d) ~= missing;
        xs = X(mask,d); % vector will all not missing values
        switch C(d)
            case 'g', mu = mean(xs);
            case 'p', mu = min(xs) - 10^-10; 
            case 'n', mu = min(xs) - 1;
            case 'c', mu = min(xs) - 1;
            case 'o', mu = min(xs) - 1;
            otherwise
                error('Unknown data type: should be g, p, c, n, c, o');
        end
        if ismember(C(d),['g','p','n'])
            suffStats{d}(1) = -mu;              % offset
            suffStats{d}(2) = 2 / max(xs-mu);   % slope
            X(mask,d) = 2 .* (xs - mu)./ max(xs-mu);
        else
            X(mask,d) = X(mask,d) - mu; % be sure that labels start at 1 (C++ implementation)
        end
    end

%     W = zeros(1,D); % vector of weights for transformation
%     for d=1:D
%         % normalize data (such that it occupies interval (0,max(X(d))]
%         if ((data.C(d) == 'n') || (data.C(d) == 'c') || (data.C(d) == 'o'))   && (min(data.X(:,d)) > 1)
%             offset = min(Xmiss(Xmiss(:,d) ~= missing,d));
%             Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 1;
%         end
%         if (data.C(d) == 'p') && (min(data.X(:,d)) > 0)
%             offset = min(data.X(:,d));
%             Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 10^-6;
%         end
%         W(d) = 2/max( Xmiss(Xmiss(:,d) ~= missing,d) );
%         
%         if ((data.C(d) == 'n' || (data.C(d) == 'c') || (data.C(d) == 'o')) && (min(data.X(:,d)) == 0))
%             Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 1;
%         elseif (data.C(d) == 'p') && (min(data.X(:,d)) == 0)
%             Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 10^-6;
%         end
%     end