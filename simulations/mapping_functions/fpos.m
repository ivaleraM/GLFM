function x = fpos(y, w)
    % Inputs:
    %   y: 
    %   w: scalar

    % TODO: make weights W input-dependent
    
    %x = log( exp( w * y ) + 1);
    
    x = 0.5 * y^2;
    
    %x = log(exp(y)+1) / w;
end
