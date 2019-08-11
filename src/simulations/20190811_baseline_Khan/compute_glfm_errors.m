function err_D = compute_glfm_errors(X_true,X_pred,X_miss,C, R)
    % miss: idxs of missings
    % X_true [N x D]
    % X_pred [N x D]
    % R [1 x D]: number of classes for each discrete dimension
    % C [1 x D]: data type per dimension
    [N,D] = size(X_true);
    if not(size(X_true) == size(X_pred))
        error('Inconsistent sizes between X_true and X_pred');
    end
    err_D = zeros(1,D);
    [n_idxs, d_idxs] = find(isnan(X_miss));
    for i=1:length(d_idxs) % for each missing value
        d = d_idxs(i);
        XT = X_true(n_idxs(i),d);
        MU = X_pred(n_idxs(i),d);
        if (C(d)=='g' )
            err_D(d)= err_D(d)+ (XT - MU)^2;
        elseif (C(d)=='p')
            err_D(d)= err_D(d)+ (XT - MU)^2;
        elseif (C(d)=='n' )
            err_D(d)= err_D(d)+ (XT - round(MU))^2;
        elseif (C(d)=='c')
            err_D(d)= err_D(d)+ (XT ~= round(MU));
        elseif (C(d)=='o')
            err_D(d)= err_D(d)+ abs(XT - round(MU))/R(d);
            %Err(p,it) =Err(p,it)+ abs(XT - MU)/range(XT(:,d));
        end
    end