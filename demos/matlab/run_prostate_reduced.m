clear
for ll=1:10
    params.simId = ll;
    params.Niter = 10000;
    params.save = 1;
    params.bias = 0;
    params.s2B = 1;
    demo_data_exploration_reduced; %_bias2;
end