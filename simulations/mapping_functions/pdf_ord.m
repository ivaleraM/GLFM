function pdf = pdf_ord(Zn,Bd,theta,s2y)
    % Function to compute pdf of an ordinal variable. It returns the whole
    % pdf, a probability vector of length R (number of categories)
    % Input parameters:
    %    Zn: [1*K], feature activation vector
    %    Bd: [K*1], feature weights (dictionary)
    % theta: [1*(R-1)]
    %   s2y: scalar, variance of pseudo-observations
    R = length(theta)+1; % number of categories
    pdf = zeros(1,R);
    for r=1:R
        if (r==1)
            a = normcdf(theta(1),Zn * Bd,s2y);
            b = 0;
        elseif (r== R)
            a = 1;
            b = normcdf(theta(r-1),Zn * Bd,s2y);
        else
            a = normcdf(theta(r),Zn * Bd,s2y);
            b = normcdf(theta(r-1),Zn * Bd,s2y);
        end
        pdf(r) = a - b;
    end
end