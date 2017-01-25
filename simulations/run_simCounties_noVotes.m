function run_simCounties_noVotes(Niter_in,Niter_out, s2Y, s2u, alpha, simId)
addpath(genpath('../Ccode/'));
%Niter_in = 1;
%Niter_out = 1;
%s2Y = 1;
%s2u = 1;
%alpha = 1;

missing=-1;

%% Selecting missing data
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

load ../databases/dataExploration/mat/counties.mat %../databases/Wine.mat

% for data exploration, remove nan
idx_unknownCounties = find( sum(isnan(data.X),2) ); 
data.X(idx_unknownCounties,:) = [];

mask_noMiss = ~isnan( data.X(:,14) );
dem_id = (data.X(mask_noMiss,14) > data.X(mask_noMiss,15)) & (data.X(mask_noMiss,14) > data.X(mask_noMiss,16));
rep_id = (data.X(mask_noMiss,15) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,15) > data.X(mask_noMiss,16));
per_id = (data.X(mask_noMiss,16) > data.X(mask_noMiss,14)) & (data.X(mask_noMiss,16) > data.X(mask_noMiss,15));

idx_to_remove = [1,3,4];
data.X(:,idx_to_remove) = []; % remove dimensions with excessive number of missings
data.C(idx_to_remove) = [];
data.cat_labels(idx_to_remove) = [];
data.ylabel(idx_to_remove) = [];


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
Zini= [dem_id, rep_id, per_id, double(rand(N,1)>0.8)];
bias = 3;
Zest = Zini';
for it=1:Niter_out
    [Zest B Theta]= IBPsampler(Xmiss,data.C,Zest',bias,W,s2Y,s2u,s2B,alpha,Niter_in,maxK,missing);
    sum(Zest')
    toc;
end
save(sprintf( '../results/counties_Niter%d_%d_s2Y%.2f_s2u%.2f_alpha%.2f.mat',Niter_in,Niter_out, s2Y, s2u, alpha) );
