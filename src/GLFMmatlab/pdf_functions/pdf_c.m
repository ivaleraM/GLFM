function pdf = pdf_c(Zp, B, params)
    %
    % Inputs:
    %   B: K*R
    %   Zp: 1*K, where K: number of latent features
    
    numMC = 1000; % number of Monte Carlo samples to approximate the Expectation
    R = size(B,2);
    pdf = zeros(1,R);
    % generate r.v's u
    u = sqrt(params.s2y) * randn(numMC,1);
    for r=1:R
        for j=1:R
            if (j==r)
                continue;
            end
            pdf(r) = pdf(r) + mean( normcdf( u + Zp*(B(r) - B(j)) ) );
        end
    end
end