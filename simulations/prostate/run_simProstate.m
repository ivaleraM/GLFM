
function run_simProstate(Niter, s2Y, s2u, alpha, simId)
       
       addpath(genpath('../src/Ccode/'));
       randn('seed',round(sum(1e5*clock)));
       rand('seed',round(sum(1e5*clock)));

       %% LOAD DATA

       load ../../datasets/mat/prostate.mat

       %% ADAPT INPUT DATA --> put bias
       drug_identifier = data.X(:,2) > 0.5;
       % remove drug levels
       data.X(:,2) = [];
       data.C(2) = [];
       data.cat_labels(2) = [];
       data.ylabel(2) = [];

       % replace missings + preprocess
       missing = -1;
       Xmiss=data.X; % Observation matrix
       Xmiss(isnan(Xmiss)) = missing;
       [N, D]= size(data.X);
       %s2Y=.5;    % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
       %s2u=.5;
       s2B=1;      % Variance of the Gaussian prior of the weigting matrices B
       %alpha=1;    % Concentration parameter of the IBP
       %Nsim=100; % Number of iterations for the gibbs sampler
       maxK= D;

       %% Inference
       tic;
       Zini= [drug_identifier, not(drug_identifier), double(rand(N,1)>0.8)];
       bias = 2;
       Zest = Zini';
       for it=1:1
           [Zest B Theta]= IBPsampler(Xmiss,data.C,Zest',bias,W,s2Y,s2u,s2B,alpha,Niter,maxK,missing);
           sum(Zest')
           toc;
       end
       save(sprintf( '../results/prostate_simId%d_%dx%d_s2Y%.2f_s2u%.2f_alpha%.2f.mat',simId,Niter_in,Niter_out, s2Y, s2u, alpha) );

end
