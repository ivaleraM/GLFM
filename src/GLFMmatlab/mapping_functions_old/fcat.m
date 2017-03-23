function x = fcat(y)
    % input argument y: [N*R]
    % output: x [N*1]
    [mm, x] = max(y,[],2);
end