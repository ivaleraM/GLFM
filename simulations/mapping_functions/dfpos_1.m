function y = dfpos_1(x, w)
    % w: scalar
        
    %y = 1./sqrt(x);
    %y = -0.5*w^0.5 * x.^(-0.5);
    y = w ./ ( 1 - exp( w .* x) );
    
end