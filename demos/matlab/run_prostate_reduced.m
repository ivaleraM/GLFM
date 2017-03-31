clear
for ll=1:5
    params.simId = ll;
    params.Niter = 5000;
    params.save = 1;
    params.s2B = 1;
    demo_data_exploration_red_bias2;
end