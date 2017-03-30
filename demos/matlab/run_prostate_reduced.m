clear
for ll=1:5
    params.simId = ll;
    params.Niter = 1000;
    params.save = 1;
    params.s2B = 1;
    demo_data_exploration_reduced;
end