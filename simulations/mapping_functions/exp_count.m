function expect = exp_count(Zn,Bd, params)
    % function to compute expectation of random variable
    
    s2y = params{1};
    fpos_1_handler = params{2};
    w = params{3}; % function to compute normalization weights w
    maxX = params{4};
    
    xin = 1:1:maxX;
    
    func = @(xin) pdf_count(xin,Zn,Bd,w, s2y, fpos_1_handler);
        
    expect = xin * func(xin)';