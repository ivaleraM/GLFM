function pdf = pdf_real(X, Zn,Bd,s2y,s2u)
    % Probability Density Function for real variables
    pdf = normpdf(X, Zn* Bd, sqrt(s2y + s2u) );