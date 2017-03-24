function [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, d)
    % Function to generate the PDF solutions corresponding to patterns in
    % Zp, and dimension d
    % Inputs:
    %   C: 1*D string with data types, D = number of dimensions
    %   Zp: P * K matrix of patterns (P is the number of patterns)
    %   B: latent feature matrix (D * K * maxR)   
    %   T: structure with parameters for mapping function
    %       T.mu: 1*D shift parameter
    %       T.w:  1*D scale parameter
    % Theta: D*maxR matrix where R is the max number of categories
    %
    % Outputs:
    %   xd: 1*numS where numS is the number of points to be considered
    %  pdf: P*numS where P is the number of patterns to consider
    
    P = size(Zp,1);
    K = size(hidden.B,2);
    if (size(Zp,2) ~= K)
        error('Incongruent sizes between Zp and hidden.B');
    end
    if (data.C == 'g') || (data.C == 'p')
        numS = 100;
        xd = linspace( min(data.X(data.X(:,d) ~= params.missing, d)), max(data.X(data.X(:,d) ~= params.missing, d)), numS);
    elseif (data.C == 'n')
        xd = min(data.X(data.X(:,d) ~= params.missing, d)):max(data.X(data.X(:,d) ~= params.missing, d));
        numS = length(xd);
    else
        xd = unique(data.X(data.X ~= missing, d));
        numS = length(xd); % number of labels for categories or ordinal data
    end
    pdf = zeros(P,numS);
    for p=1:P
        switch data.C(d)
            case 'g', pdf(p,:) = pdf_g(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.w(d), params);
            case 'p', pdf(p,:) = pdf_p(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.w(d), params);
            case 'n', pdf(p,:) = pdf_n(x,Zp, squeeze(hidden.B(d,:,1))', params);
            case 'c', pdf(p,:) = pdf_c(Zp, squeeze(hidden.B(d,:,1:hidden.R(d))), params);
            case 'o', pdf(p,:) = pdf_o(Zp, squeeze(hidden.B(d,:,1))', hidden.Theta(d,:), params);
            otherwise
                error('Unknown data type');
        end
        if (sum(isnan(X_map(:,d))) > 0)
            error('Some values are nan!');
        end
    end
end