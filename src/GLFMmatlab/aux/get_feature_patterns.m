function [patterns, C, L] = get_feature_patterns(Z,varargin)
% function to return matrix with all patterns dim = (numPatterns*numFeatures)
% Inputs:
%   - Z: N*K matrix
% Outputs:
%   - patterns: (numPatterns*numFeatures) matrix
%   - C is an assignment vector, with the pattern id for each observation
%   - L is a numPatterns*1 vector with the number of samples per pattern

% average
N = size(Z,1);
C = zeros(N,1);

if isempty(varargin)
    patterns = unique(Z,'rows'); % all existing patterns
else
    patterns = varargin{1};
end

L = zeros(1,size(patterns,1));
for r=1:size(patterns,1) % for each pattern
    mask = sum( ~xor(repmat(patterns(r,:),size(Z,1),1), Z),2) == size(Z,2); % which patients have that pattern
    C(mask) = r;
    L(r) = sum(mask);
    fprintf('%d. %s: %d\n', r, num2str(patterns(r,:)), L(r));
end
