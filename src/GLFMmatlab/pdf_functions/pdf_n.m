function pdf = pdf_n(x,Zp, B, params)
    %
    % Inputs:
    %   B: K*R
    %   Zp: 1*K, where K: number of latent features
    
    pdf = cdf( (f_n_1(x+1,mu,w) - Zp*B )/sqrt(params.s2y) ) - ...
        cdf( (f_n_1(x,mu,w) - Zp*B )/sqrt(params.s2y) );
end