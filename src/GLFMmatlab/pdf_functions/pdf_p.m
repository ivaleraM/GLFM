function pdf = pdf_p(x,Zp, B, mu, w, s2Y, params)

% IMPORTANT: func = 2 does not have f_p_1 function implemented
% func 2 = deprecated

    pdf = 1./(2*pi*sqrt(s2Y + params.s2u)) * ...
        exp( -1/(2*(s2Y + params.s2u)) .* ...
        (f_p_1(x, mu, w) - Zp * B).^2 ) .* abs(df_p_1(x, mu, w));
end