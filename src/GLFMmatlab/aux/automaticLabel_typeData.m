function S = automaticLabel_typeData(X)
% possible data types:
%   g - gaussian distributed (in R)
%   p - positive gaussian (in R+)
%   n - integer (in N)
%   c - categorical (in category set F)
%   o - ordinal (in ordered category set O)

S = cell(1,size(X,2));

% detect g
maskG = logical(sum(X < 0,1));
S(maskG) = {'g'};

% detect p
maskP = not(sum(X<0,1)) & sum(mod(X,1)~=0); % column does not have neg numbers & has some floating contribution
S(maskP) = {'p'};

% detect c (all categoricals in our case are binary, because we have
% expanded them)
maskC = not( sum((X ~= 0) & (X~= 1),1) ); % column has only ones and zeros
S(maskC) = {'c'};

% all the rest, count data
already = maskG + maskP + maskC;
if sum(already > 1)
    error('Should never overlap')
end
idxs = setdiff(1:size(X,2),find(already));
S(idxs) = {'n'};

S = cell2mat(S);