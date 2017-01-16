function pdf = pdf_pos(X,Z,B,W,s2y)

    pdf = 1/sqrt(2*pi*(1+s2y)) .* exp( -1./(2.*(s2y+1)) .* (fpos_1(X,W) - Z * B).^2 ) ...
        .* abs( dfpos_1(X) );
