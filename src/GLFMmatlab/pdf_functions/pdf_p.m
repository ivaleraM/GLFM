function pdf = pdf_p(x,Zp, Bd, mu, w, s2Y, params)
    %
    % Likelihood function for positive data
    % Eq. (2) in the paper:
    % "General Latent Feature Models for Heterogeneous Datasets"
    %
    % Inputs:
    %     x: positive observation
    %    Bd: K*1  weight vector for individual dimension d
    %    Zp: 1*K, binary vector of feature assignment,
    %          where K: number of latent features
    %   s2Y: Gaussian noise variance
    %  mu,w: Hyper-parameters linked to the transformation described in the paper

    pdf = 1./sqrt(2*pi*(s2Y + params.s2u)) * ...
        exp( -1/(2*(s2Y + params.s2u)) .* ...
        (f_p_1(x, mu, w) - Zp * Bd).^2 ) .* abs(df_p_1(x, mu, w));
end
