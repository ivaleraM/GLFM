function pdf = pdf_pos(X,Zn,Bd,w,s2y,s2u)
    % Probability density function for positive real variables
    % Inputs:
    %   X: values in which to evaluate pdf
    %  Zn: [1*K] vector
    %  Bd: [K*1] vector
    %   w: scalar: normalization weight % TODO: Review how to automatize
    % s2y: scalar, variance of pseudo-observations
    % s2u: scalar, variance of auxiliary noise
    
    pdf = 1/sqrt(2*pi*(s2u+s2y)) .* exp( -1/(2*(s2y+s2u)) * (fpos_1(X,w) - Zn*Bd).^2 ) ...
        .* abs( dfpos_1(X,w) );
