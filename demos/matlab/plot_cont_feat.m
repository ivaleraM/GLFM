function V = plot_cont_feat(X, labels, varargin)
% X is 2*D array with sctive vector + ref. vector

V = log2( X(1,:) ./ X(2,:) );
figure; plot(V, 'Linewidth', 2);
xticks(1:size(X,2));
xticklabels(labels);
xtickangle(45);
grid;

if ~isempty(varargin)
    title(varargin{1});
end