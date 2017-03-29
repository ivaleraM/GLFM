function pdf = pdf_g(x,Zp, B, w, s2Y, params)
    pdf = normpdf(x, Zp * B, sqrt(s2Y + params.s2u)) .* w;
end