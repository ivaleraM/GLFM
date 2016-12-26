% MAIN TOY SCRIPT (to compare behavior of Matlab and Python)
addpath('../Ccode/');

X = [1.0 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15];
Z = [1.0, 0; 1 1; 1 1];
C = 'PPPPP';

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

disp('SUCCESSFUL');


