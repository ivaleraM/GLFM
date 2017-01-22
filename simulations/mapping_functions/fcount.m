function x = fcount(y,w)
    % Inputs:
    %   x: [N*D]
    
    x = floor( fpos(y, w) );
end