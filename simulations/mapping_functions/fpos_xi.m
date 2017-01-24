function x = fpos_xi(y, w)
    % Inputs:
    %   y: 
    %   w: scalar

    % TODO: make weights W input-dependent
    
    %x = log( exp( w * y ) + 1);
    
    x =  y^2 /w;
    
    %x = log(exp(y)+1) / w;
end
