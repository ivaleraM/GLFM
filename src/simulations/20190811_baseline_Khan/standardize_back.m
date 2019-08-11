function X_stand = standardize_back(X_pred, data)
    X_stand.continuous = X_pred.continuous;
    X_stand.discrete = X_pred.discrete;

    if (~isfield(data,'mu') || ~isfield(data,'std_dev')) || ~isfield(data,'cell_nCat')
        error('Standardization info missing');
    end
    % revert standardization continuous variables
    D = size(X_pred.continuous,1);    
    for d=1:D
        V = X_pred.continuous(d,:);
        V(not(isnan(V))) = V(not(isnan(V))) * data.std_dev(d) + data.mu(d);
        X_stand.continuous(d,:) = V;
    end
    
    % revert standardization discrete variables
    Dcat = size(X_pred.discrete,1);
    for d=1:Dcat
        V = X_pred.discrete(d,:);
        for r=1:length(data.cell_nCat{d})
            cat = data.cell_nCat{d}(r);
            V(X_pred.discrete(d,:) == r) = cat;
        end
        X_stand.discrete(d,:) = V;
    end