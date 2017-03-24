function pdf = pdf_o(Zp, B, theta, params)
    % Inputs:
    %   theta: 1*R
    %   Zp: 1*K
    %   B: K*D
    R = length(theta);
    pdf = zeros(1,R);
    for r=1:R
        if (r==1)
            a = normcdf(theta(r), Zp*B, sqrt(params.s2Y) );
            b = 0;
        elseif (R==R)
            a = normcdf(theta(r), Zp*B, sqrt(params.s2Y) );
            b = normcdf(theta(r-1), Zp*B, sqrt(params.s2Y) );
        end
        pdf =  a - b;
    end
end