function res = max_not_missing(V, missing_val)

    res = zeros(1,size(V,2));
    for d=1:size(V,2)
        res(d) = max(V( V(:,d)~=missing_val,d));
    end