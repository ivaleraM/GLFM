% MAIN TOY SCRIPT (to compare behavior of Matlab and Python)
addpath(genpath('../Ccode/'));
clear

X = [1.0 1 -0.3 1; 6 2 3.8 23; 11 3 4.1 4];
Z = [1.0, 0; 1 1; 1 1];
C = 'PoGN';

disp('First, in Python')
X
Z

%print '\nX flags'
%print X.flags
%print '\n'
%print X.transpose().flags
%
%print '\nZ flags'
%print Z.flags
%print '\n'
%print Z.transpose().flags

% Define default values
bias = 1; s2Y = 1.0; s2B = 1.0; alpha = 1.0; Nsim = 50; maxK = 10; missing = -1;

disp('Now, inside C');
%wrapper_IBPsampler(X,C,Z)
[Zest B Theta] = IBPsampler(X,C,Z,bias, s2Y, s2B, alpha, Nsim, maxK, missing)
disp('Returned from C');

disp('SUCCESSFUL');


