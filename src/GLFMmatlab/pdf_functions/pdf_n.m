function pdf = pdf_n(x,Zp, Bd, mu, w, s2Y, params)
    %
    % Likelihood function for count data
    % Eq. (8) in the paper:
    % "General Latent Feature Models for Heterogeneous Datasets"
    %
    % Inputs:
    %    x: count-data observation
    %   Bd: K*1  weight vector for a particular dimension d
    %   Zp: 1*K, binary vector of feature assignment,
    %       where K: number of latent features
    %  s2Y: Gaussian noise variance  for u
    % mu,w: Hyper-parameters linked to the transformation described in the paper

    pdf = normcdf(f_n_1(x+1,mu,w), Zp*Bd,sqrt(s2Y) ) - ...
        normcdf(f_n_1(x,mu,w), Zp*Bd, sqrt(s2Y) );
end
