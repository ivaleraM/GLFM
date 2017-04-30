function x = f_c(y)
    % transformation function for categorical data
    % (optional) eventual input = label
    % Input y: N*R where N is the number of observations and R num of
    % categories
    [a,x] = max(y,[],2);
end
