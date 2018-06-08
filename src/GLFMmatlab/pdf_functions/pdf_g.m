function pdf = pdf_g(x,Zp, B, mu, w, s2Y, params)

    % Likelihood function for real-valued data
    % Eq. (2) in the paper:
    % "General Latent Feature Models for Heterogeneous Datasets"
    %
    %    x: real observation
    %    B: K*1  weight vector for feature k
    %   Zp: 1*K, binary vector of feature assignment,
    %       where K: number of latent features
    %  s2Y: Gaussian noise variance for pseudo-observations
    %           params.s2u: variance for additive Gaussian noise
    %           in the transformation x = f(y+u)
    %    w: Hyper-parameter linked to the transformation described in the paper


     df_1 = @(x) w*(x-mu);
    pdf = normpdf( df_1(x) , Zp * B, sqrt(s2Y + params.s2u)) .* w;
end
