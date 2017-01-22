function x = freal(y, s2u)
    % Mapping function for real values
    x = y + sqrt(s2u) * randn(size(y));
end