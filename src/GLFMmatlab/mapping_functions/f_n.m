function x = f_n(y, mu, w)
    % transformation function for count data
    % Y -> X (from pseudo-obversations to data)
    if (w == 0)
        error('scaling factor should never be 0');
    end
    x = floor( f_p(y, mu, w) );
    %if sum(x == 0)>0 % TODO: verify if that would be a problem % DONE, OK
    %    warning('when transforming y to x, x has been rounded to zero (count data)');
    %end
end