function res = std_not_isnan(V,dim)

if (dim == 1)
    V = V';
end

res = zeros(1,size(V,2));
for d=1:size(V,2)
    res(d) = std(V(~isnan(V(:,d)),d));
end