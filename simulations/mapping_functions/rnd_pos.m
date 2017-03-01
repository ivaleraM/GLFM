function samples = rnd_pos(Zn,Bd,numSamples,params)
    % function to get random samples from the distribution

    s2y = params{1};
    s2u = params{2};
    fpos_1_handler = params{3};
    dfpos_1_handler = params{4};
    w = params{5}; % function to compute normalization weights w
    maxX = params{6};

    func = @(xin) log( pdf_pos(xin,Zn,Bd,w, s2y ,s2u, fpos_1_handler, dfpos_1_handler) );
    a = 10^-6;
    b = maxX;
    domain = [a,b+5];

    samples = ars(func, a, b, domain, numSamples, []);
    samples = exp(samples);
