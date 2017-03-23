function res = min_not_missing(V, missing_val)

    res = zeros(1,size(V,2));
    for d=1:size(V,2)
        res(d) = min(V( V(:,d)~=missing_val,d));
    end
