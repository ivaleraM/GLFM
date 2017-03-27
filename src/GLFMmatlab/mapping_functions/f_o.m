function x = f_o(y, Theta)
    % transformation function for ordinal data
    x = zeros(size(y)); % column vector
    for j=1:length(Theta)
        if (j == 1)
            mask = (y <= Theta(1));
        else
            mask = (y > Theta(j-1)) & (y <= Theta(j));
        end
        x( mask ) = j;
    end
    x ( x == 0 ) = length(Theta); % last ordinal category
end