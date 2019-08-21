function y = f_p_1(x, mu, w)
    % transformation function for positive data
    % X -> Y (from data to pseudo-obversations)
    if (w == 0)
        error('scaling factor should never be 0');
    end
    %%y = log( exp(w.*(x-mu) - 1) ); previous wrong version
    y = log( exp(w.*(x-mu) )- 1 );
end