function T = compute_jaccard_simMatrix(Z)

K = size(Z,2);
N = size(Z,1);

T = zeros(K);
for k=1:K
    for j=k:K
        T(k,j) = sum(Z(:,k) & Z(:,j))./ sum( Z(:,k) | Z(:,j) );
        T(j,k) = T(k,j);
    end
end