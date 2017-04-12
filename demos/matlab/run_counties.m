for ll=1:5
    params.simId = ll;
    params.bias = 1;
    params.Niter = 5000;
    params.save = 1;
    demo_counties;
    
    params.bias = 0;
    demo_counties;
end