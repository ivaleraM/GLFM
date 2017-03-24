function pdf = pdf_g(x,Zp, B, w, params)
    pdf = normpdf(x, Zp * B, sqrt(params.s2Y + params.s2u)) .* w;
end