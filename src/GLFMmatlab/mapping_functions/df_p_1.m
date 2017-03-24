function y = df_p_1(x, mu, w)
    % derivative of transformation function for positive data
    % X -> Y (from data to pseudo-obversations)
    if (w == 0)
        error('scaling factor should never be 0');
    end
    y = exp(w.*x) ./ ( exp(w.*(x-mu) - 1) );