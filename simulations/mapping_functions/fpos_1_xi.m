function y = fpos_1_xi(x,w)
    % w: scalar
    
    y = sqrt(w*x);
    % y = log( exp( x ) + 1) / w;
    % y = log( exp( w*x ) -1 );
end