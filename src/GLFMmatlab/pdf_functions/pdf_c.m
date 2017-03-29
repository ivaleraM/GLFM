function pdf = pdf_c(Zp, B, s2Y)
    %
    % Inputs:
    %   B: K*R
    %   Zp: 1*K, where K: number of latent features
    
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