function y = dfpos_1(x)
    
    W = repmat( 2 ./ max(x,[],1) , size(x,1), 1);
    y = W ./ ( 1 - exp( W .* x) );
    
end