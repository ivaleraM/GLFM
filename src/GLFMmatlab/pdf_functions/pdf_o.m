function pdf = pdf_o(Zp, Bd, theta, s2Y)
    %
    % Likelihood function for ordinal data
    % Eq. (6) in the paper:
    % "General Latent Feature Models for Heterogeneous Datasets"
    %
    % Inputs:
    %    Bd: K*1  weight vector for dimension d
    %    Zp: 1*K, binary vector of feature assignment,
    %          where K: number of latent features
    %   s2Y: Gaussian noise variance  for u
    % theta: 1xR auxiliary thresholds, see description in the paper

    R = length(theta)+1;
    pdf = zeros(1,R);
    for r=1:R
        if (r==1)
            a = normcdf(theta(r), Zp*Bd, sqrt(s2Y) );
            b = 0;
        elseif (r==R)
            a = 1;
            b = normcdf(theta(r-1), Zp*Bd, sqrt(s2Y) );
        else
            a = normcdf(theta(r), Zp*Bd, sqrt(s2Y) );
            b = normcdf(theta(r-1), Zp*Bd, sqrt(s2Y) );
        end
        pdf(r) =  a - b;
    end
end
