function res = mean_not_isinf(V,dim)

if (dim == 0)
    res = zeros(1,size(V,2));
    for d=1:size(V,2)
        res(d) = mean(V(~isinf(V(:,d)) & ~isnan(V(:,d)) ,d));
    end
else
    error('Implement!')
end