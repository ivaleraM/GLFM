function pdf = pdf_o(Zp, B, theta, params)
    % Inputs:
    %   theta: 1*R
    %   Zp: 1*K
    %   B: K*D
    R = length(theta)+1;
    pdf = zeros(1,R);
    for r=1:R
        if (r==1)
            a = normcdf(theta(r), Zp*B, sqrt(params.s2Y) );
            b = 0;
        elseif (r==R)
            a = 1;
            b = normcdf(theta(r-1), Zp*B, sqrt(params.s2Y) );
        else
            a = normcdf(theta(r), Zp*B, sqrt(params.s2Y) );
            b = normcdf(theta(r-1), Zp*B, sqrt(params.s2Y) );
        end
        pdf(r) =  a - b;
    end
end