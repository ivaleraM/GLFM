function data_transformed = transform_dataset(data,params)

    % replace missings
    data_transformed = data;
    data_transformed.X(isnan(data_transformed.X)) = params.missing;

    % change labels for categorical and ordinal vars such that > 0
    %V_offset = zeros(1,D);
    [N,D] = size(data.X);
    for d=1:D
        if (data_transformed.C(d) == 'c') || (data_transformed.C(d) == 'o')
            mask = data_transformed.X(:,d) ~= params.missing;
            %V_offset(d) = min( data_transformed.X(mask,d) );
            %data_transformed.X(mask,d) = data_transformed.X(mask,d) - V_offset(d) + 1;
            uniqueVal= unique(data_transformed.X(mask,d));
            Xaux=[];
            for i=1:length(uniqueVal)
                Xaux(data_transformed.X(:,d)==uniqueVal(i)) = i;
            end
            Xaux(~mask)=params.missing;
            data_transformed.X(:,d)= Xaux;
        end
    end

    % eventually, apply external transform
    for r=1:size(data_transformed.X,2)
        if ~isempty(params.t{r})
            data_transformed.X(:,r) = params.t_1{r}(data_transformed.X(:,r)); % work in logarithm space better
            data_transformed.C(r) = params.ext_dataType{r};
        end
    end