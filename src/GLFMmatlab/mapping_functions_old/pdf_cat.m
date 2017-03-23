function pdf = pdf_cat(Zn,Bd,s2u,R)
    % Function to compute pdf of an ordinal variable. It returns the whole
    % pdf, a probability vector of length R (number of categories)
    % Input parameters:
    %    Zn: [1*K], feature activation vector
    %    Bd: [K*1], feature weights (dictionary)
    %   s2u: scalar, variance of auxiliary noise
    %     R: number of categories

    pdf = zeros(1,R);
    numMC_samples = 100;
    % TODO: perform pdf computation in log space ?
    for r=1:R
        tmp = ones(1,numMC_samples);
        % we compute the expectation using Monte Carlo samples
        uV = sqrt(s2u) * randn(1,numMC_samples); % mean = 0 % TODO: check that u does not depend on j
        for j=1:R
            if (j==r)
                continue;
            end
            tmp = tmp .* normcdf(uV + repmat( Zn * (Bd(:,r)-Bd(:,j)), 1, numMC_samples) );
        end
        pdf(r) = mean(tmp,2);
    end
