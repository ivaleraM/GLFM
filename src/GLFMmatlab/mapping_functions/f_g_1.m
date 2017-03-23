function y = f_g_1(x, mu, w)
    % transformation function for real-valued data
    % X -> Y (from data to pseudo-obversations)
    if (w == 0)
        error('scaling factor should never be 0');
    end
    y = w .* (x - mu);
end