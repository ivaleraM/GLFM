function data = standardize(data)
    
    if ( sum(isnan(data.continuous(:))) || sum(isnan(data.discrete(:))) )
        error('Standardization should be done before applying any missing to avoid diff cat in test/train')
    end
    % standardize continuous variables
    D = size(data.continuous,1);
    data.mu = zeros(1,D);
    data.std_dev = zeros(1,D);
    for d=1:D
        data.mu(d) = mean_not_isnan(data.continuous(d,:),1);
        data.std_dev(d) = std_not_isnan(data.continuous(d,:),1);
        if isnan(data.std_dev(d)) | isnan(data.mu(d)) | isinf(data.std_dev(d)) | isinf(data.mu(d))
            error('Infinite or unknowm mu,std: standardization of data not possible')
        end
        V = data.continuous(d,:);
        V(not(isnan(V))) = (V(not(isnan(V))) - data.mu(d))*1.0/data.std_dev(d);
        data.continuous(d,:) = V;
    end
    
    % standardize discrete variables
    [Dcat,N] = size(data.discrete);
    data.cell_nCat = cell(1,Dcat);
    for d=1:Dcat
        data.cell_nCat{d} = unique(data.discrete(d,:));
        data.cell_nCat{d} = data.cell_nCat{d}(not(isnan(data.cell_nCat{d})));
        V = data.discrete(d,:);
        for r=1:length(data.cell_nCat{d})
            cat = data.cell_nCat{d}(r);
            V(data.discrete(d,:) == cat) = r;
        end
        data.discrete(d,:) = V;
    end
    data.X = data.X_true;
    data.X(:,data.idxs_continuous) = data.continuous';
    data.X(:,data.idxs_discrete) = data.discrete';