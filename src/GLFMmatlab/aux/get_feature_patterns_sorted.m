function [patterns, C, L] = get_feature_patterns_sorted(Z,varargin)
% function to return matrix with all patterns dim = (numPatterns*numFeatures)
% Inputs:
%   - Z: N*K matrix
%   - (optional 1) fid: file identifier where to print description (default: print
%   in command window with fid=1)
%   - (optional 2) patterns: P*K (complete list of patterns, this might be
%   useful if we already know the list of patterns and Z is a subset of the
%   matrix, so not all patterns might be present in that subset).
% Outputs:
%   - patterns: (numPatterns*numFeatures) matrix
%   - C is an assignment vector, with the pattern id for each patient (numPatients*1)

% average
N = size(Z,1);
D = size(Z,2); % number of features
C = zeros(N,1);

if isempty(varargin)
    fid = 1;
    patterns = unique(Z,'rows'); % all existing patterns
elseif length(varargin) == 1
    fid = varargin{1};
    patterns = unique(Z,'rows');
else
    fid = varargin{1};
    patterns = varargin{2};
end

L = zeros(1,size(patterns,1));
for r=1:size(patterns,1) % for each pattern
    mask = sum( ~xor(repmat(patterns(r,:),size(Z,1),1), Z),2) == size(Z,2); % which patients have that pattern
    C(mask) = r;
    L(r) = sum(mask);
end

% sort patterns
[L, I] = sort(L, 'descend');
patterns = patterns(I,:);

L = zeros(1,size(patterns,1));
for r=1:size(patterns,1) % for each pattern
    mask = sum( ~xor(repmat(patterns(r,:),size(Z,1),1), Z),2) == size(Z,2); % which patients have that pattern
    C(mask) = r;
    L(r) = sum(mask);
    fprintf(fid, '%d. %s: %d\n', r, num2str(patterns(r,:)), L(r));    
end
