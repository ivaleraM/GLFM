%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Description of variables and structs for IBPsampler %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
% X:       Observation matrix of size $N x D$ where N is number of observations and
%          D is the dimension of the observations. Here missing data should
%          be introduced as -1, and categorical and ordinal data should take
%          values in {1,2, ..., R}.
% C:       Char array of length D, indicating the type of variable in each
%          column of the observation matrix X. Here 'g' indicates Real
%          variable, 'p' positive real variable, 'n' count data, 'o' ordinal
%          data and 'c' categorical data. % TODO: Complete dimensions!
% bias:    Number of colums of Z not to be sampled
% Zini:    initialization of the IBP matrix. % TODO: Complete dimensions!
% s2Y:     Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
% s2B:     Variance of the Gaussian prior of the weigting matrices B
% alpha:   Concentration parameter of the IBP
% Nsim:    Number of iterations for the gibbs sampler
% maxK:    Maximum number of latent fatures allowed
% missing: Value to which each missing data is encoded to

%% Outputs
% Zest:  Inferred IBP matrix
% B:     Weighting matrix. Size [D x Kest x maxR]
% Theta: Inferred thresholds for ordinal data. Size [D x maxR].

%% Compute missing data
% Apply the corresponding f(Y) -- See function matrix_completion.m
