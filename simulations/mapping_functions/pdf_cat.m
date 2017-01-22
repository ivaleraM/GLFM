function pdf = pdf_cat(Zn,B,s2u)
    % Function to compute pdf of an ordinal variable. It returns the whole
    % pdf, a probability vector of length R (number of categories)
    % Input parameters:
    %    Zn: [1*K], feature activation vector
    %    Bd: [K*1], feature weights (dictionary)
    %   s2u: scalar, variance of auxiliary noise

    pdf = zeros(1,maxR);
    numMC_samples = 100;
    % TODO: perform pdf computation in log space ?
    for r=1:maxR
        tmp = ones(1,numMC_samples);
        % we compute the expectation using Monte Carlo samples
        uV = sqrt(s2u) * randn(1,numMC_samples); % mean = 0 % TODO: check that u does not depend on j
        for j=1:maxR
            if (j==r)
                continue;
            end
            tmp = tmp .* normcdf(uV + repmat( Zn * (B(:,d,r)-B(:,d,j)), 1, numMC_samples) );
        end
        pdf(r) = mean(tmp,2);
    end
