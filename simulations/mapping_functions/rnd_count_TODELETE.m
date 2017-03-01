function x = rnd_count_TODELETE(Zn,Bd,numSamples, params)
    % function to get random samples from the distribution

    s2y = params{1};
    fpos_1_handler = params{2};
    w = params{3}; % function to compute normalization weights w
    maxX = params{4};

    func = @(xin) log( pdf_count(xin,Zn,Bd,w, s2y, fpos_1_handler) );

    a = 1;
    b = maxX;
    domain = [a,b+5];

    samples = ars(func, a, b, domain, numSamples, []);
    samples = round(exp(samples));    

    %[mu_prev, Puntos] = ars(@lc_lpdf, a-10, b, domain, 1, Puntos, varargin);
