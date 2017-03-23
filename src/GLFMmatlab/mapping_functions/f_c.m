function y = f_c(x)
    % transformation function for categorical data
    % (optional) eventual input = label
    % Input x: P*R where P is the number of observations and R num of
    % categories
    [a,y] = max(x,[],2);
end