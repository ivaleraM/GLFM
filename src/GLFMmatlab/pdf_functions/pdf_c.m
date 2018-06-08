function pdf = pdf_c(Zp, B, s2Y)
    %
    % Likelihood function for categorical data
    % Eq. (4) in the paper:
    % "General Latent Feature Models for Heterogeneous Datasets"
    %
    % Inputs:
    %   B: K*R   weight matrix, where R: number of categories
    %   Zp: 1*K, binary vector of feature assignment,
    %       where K: number of latent features
    %   s2Y: Gaussian noise variance  for u


    numMC = 1000; % number of Monte Carlo samples to approximate the Expectation
    R = size(B,2);
    pdf = zeros(1,R);
    % generate r.v's u
    u = sqrt(s2Y) * randn(numMC,1);
    for r=1:R
        aux = ones(numMC,1);
        for j=1:R
            if (j==r)
                continue;
            end
            aux = aux .* normcdf( u + Zp*(B(:,r) - B(:,j)) );
        end
        pdf(r) = mean(aux);
    end
    pdf = pdf ./ sum(pdf);
end
