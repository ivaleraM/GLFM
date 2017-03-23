function x = f_n(y, mu, w)
    % transformation function for count data
    % Y -> X (from pseudo-obversations to data)
    if (w == 0)
        error('scaling factor should never be 0');
    end
    x = floor( log( exp(y) + 1 )./w + mu );
    if sum(x == 0)>0 % TODO: verify if that would be a problem
        error('when transforming y to x, x has been rounded to zero (count data)');
    end
end