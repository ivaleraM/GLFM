function [numPatients, mean_survival, median_survival] = get_patterns_info(Y,patterns,C)
    % Given the matrix of patterns, it returns 3 vectors with number of
    % patients, mean and median survival per pattern.
    % C is an assignment vector, with the pattern id for each patient

    %vals = unique(C);
    P = size(patterns,1); %length(vals);
%     if (P ~= size(patterns,1))
%         error('P and patterns are inconsistent: Should have equal size');
%     end
    numPatients = nan*ones(1,P);
    mean_survival = nan*ones(1,P);
    median_survival = nan*ones(1,P);

    for r=1:P
        mask = find(C == r);
        numPatients(r) = length(mask);
        Nboost = numPatients(r); %floor(numPatients(r) * 0.7);
        if (Nboost == 0)
            continue;
        end
        Ym = Y(mask(randi(numPatients(r),Nboost,10000) ) ); % 1 time bootstrap
        Ymean = mean(Ym); %mean_not_isnan(Ym,0);
        Ymedian = median(Ym); %median_not_isnan(Ym,0);
        mean_survival(r) = mean_not_isnan(Ymean',0);
        median_survival(r) = mean_not_isnan(Ymedian',0);
    end
    