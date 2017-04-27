for ll=4:5
    params.simId = ll;
    params.bias = 0;
    params.Niter = 10000;
    params.save = 1;
    demo_data_exploration_counties;
    
    %params.bias = 1;
    %demo_data_exploration_counties;
end