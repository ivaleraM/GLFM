function res = median_not_isinf(V,dim)

if (dim == 0)
    res = zeros(1,size(V,2));
    for d=1:size(V,2)
        res(d) = median(V(~isinf(V(:,d)),d));
    end
else
    error('Implement!')
end