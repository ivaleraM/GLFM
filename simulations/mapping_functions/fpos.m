function [y,W] = fpos(x)

    W = repmat( 2 ./ max(x,[],1) , size(x,1), 1);
    y = log( exp( W .* x ) + 1);
end