function pdf = pdf_g(x,Zp, B, mu, w, s2Y, params)


%     %mask = (x ~= params.missing);
%     %x = (x - mean(x(mask))) ./std(x(mask));    
     df_1 = @(x) w*(x-mu);
%    if (w == 0)
%        error('weight w should never be equal to zero');
%    end
%    x = x/w + mu;

    pdf = normpdf( df_1(x) , Zp * B, sqrt(s2Y + params.s2u)) .* w;
end
