function [patterns, C] = get_feature_patterns(Z,varargin)
% function to return matrix with all patterns (numPatterns*numFeatures)
% C is an assignment vector, with the pattern id for each patient
% (numPatients*1)

% average
N = size(Z,1);
D = size(Z,2); % number of features
C = zeros(N,1);

if isempty(varargin)
    patterns = unique(Z,'rows'); % all existing patterns
    %patterns = de2bi(0:(D-1),'left-msb');
else
    patterns = varargin{1};
end


for r=1:size(patterns,1) % for each pattern
    mask = sum( ~xor(repmat(patterns(r,:),size(Z,1),1), Z),2) == size(Z,2); % which patients have that pattern
    C(mask) = r;
end