function y = dfpos_1(x, w)
    % w: scalar
    
    %W = repmat( 2 ./ max(x,[],1) , size(x,1), 1);
    y = w ./ ( 1 - exp( w .* x) );
    
end