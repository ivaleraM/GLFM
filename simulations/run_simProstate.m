
function run_simProstate(Niter_in,Niter_out, s2Y, s2u, alpha, simId)

        addpath(genpath('../Ccode/'));

       missing=-1;

       %% Selecting missing data
       randn('seed',round(sum(1e5*clock)));
       rand('seed',round(sum(1e5*clock)));

          load /Users/melanie/Desktop/GLFM/databases/dataExploration/mat
       %load ../databases/dataExploration/mat/prostate.mat %../databases/Wine.mat

       drug_identifier = data.X(:,2) > 0.5;

       % remove drug levels
       data.X(:,2) = [];
       data.C(2) = [];
       data.cat_labels(2) = [];
       data.ylabel(2) = [];

       N = size(data.X,1);
       D = size(data.X,2);

       Xmiss=data.X;        % Observation matrix

       [N, D]= size(data.X);
       %s2Y=.5;    % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
       %s2u=.5;
       s2B=1;      % Variance of the Gaussian prior of the weigting matrices B
       %alpha=1;    % Concentration parameter of the IBP
       %Nsim=100; % Number of iterations for the gibbs sampler
       maxK= D;

       Xmiss(isnan(Xmiss)) = missing;
       W = zeros(1,D); % vector of weights for transformation
       for d=1:D
           % normalize data (such that it occupies interval (0,max(X(d))]
           if ((data.C(d) == 'n') || (data.C(d) == 'c') || (data.C(d) == 'o'))   && (min(data.X(:,d)) > 1)
               offset = min(Xmiss(Xmiss(:,d) ~= missing,d));
               Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 1;
           end
           if (data.C(d) == 'p') && (min(data.X(:,d)) > 0)
               offset = min(data.X(:,d));
               Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 10^-6;
           end
           W(d) = 2/max( Xmiss(Xmiss(:,d) ~= missing,d) );

           if ((data.C(d) == 'n' || (data.C(d) == 'c') || (data.C(d) == 'o')) && (min(data.X(:,d)) == 0))
               Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 1;
           elseif (data.C(d) == 'p') && (min(data.X(:,d)) == 0)
               Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 10^-6;
           end

       end


       %% Inference
       tic;
       Zini= [drug_identifier, not(drug_identifier), double(rand(N,1)>0.8)];
       bias = 2;
       Zest = Zini';
       for it=1:Niter_out
           [Zest B Theta]= IBPsampler(Xmiss,data.C,Zest',bias,W,s2Y,s2u,s2B,alpha,Niter_in,maxK,missing);
           sum(Zest')
           toc;
       end
       save(sprintf( '../results/prostate_simId%d_%dx%d_s2Y%.2f_s2u%.2f_alpha%.2f.mat',simId,Niter_in,Niter_out, s2Y, s2u, alpha) );

end
