%% Data Exploratory Anaylisis for Prostate DB
addpath(genpath('./aux/'));
load('../results/tmp_prostate.mat');
data.cat_labels{1} = {'3', '4'};

data.ylabel_long = {'Stage', 'Drug level', 'Months of Follow-up', 'Status', 'Age in years', ...
    'Weight Index', 'Type of activity', 'History of Cardiovascular Disease', ...
    'Systolic Blood Pressure/10', 'Diastolic Blood Pressure/10', ...
    'Electrocardiogram', 'Serum Hemoglobin (g/100ml)', 'Size of Primary Tumor (cm^2)', ...
    'Combined Index of Stage and Hist. Grade', 'Serum Prostatic Acid Phosphatase', ...
    'Bone Metastases'};

Z = Zest';
%Z(:,2) = not(Z(:,2);
drug_identifier = data.X(:,2) > 0.5;

%[patterns, C] = get_feature_patterns([drug_identifier, Z]);
%print_patterns(1,data.X,patterns,C);
%Z(:,2) = [];
[patterns, C] = get_feature_patterns(Z);
print_patterns(1,data.X(:,13),patterns,C,data.ylabel_long{13});

%G1 = (C < 5);
%G2 = (C == 5) | (C == 6) | (C == 7);

maskD = (data.C ~= 'c' & data.C ~= 'o'); % all dimensions excluding categorical and ordinal
figure(1); subplot(3,1,1:2); hold off
colors = {'k--', 'b', 'r', 'c', 'm'};
for f=2:size(Z,2) % for each feature
    %fprintf('Feature %d\n', f);
    x1 = data.X(Z(:,f) == 1,maskD);
    x2 = data.X(Z(:,f) == 0,maskD);
    V = log2( mean_not_isnan(x1,0) ./ (mean_not_isnan(x2,0) + 10^-6) );
    z = 1:sum(maskD);
    %SS = std_not_isinf( tmp,0); % FLAG, DIM
    MM = V;%mean_not_isinf(V,0);
    MM(MM == -Inf) = -6;
    MM(MM == +Inf) = +6;
    %f = [MM+2*SS; flipdim(MM-2*SS,1)];
    h1 = plot(z, MM, colors{f}, 'Linewidth',2); hold on
    %legend(sprintf('mean value of log_2(%s/%s)',settings.groups.g_star_label,settings.groups.g_b_label));
    %h = fill([z; flipdim(z,1)], f, 'k'); %[7 7 7]/8)
    %set(h,'facealpha',0.2)
    
    %pause;
end
set(gca,'XTick',[])
set(gca,'XTicklabel',[])
set( gca, 'XTick', 1:sum(maskD))
set( gca,'XTickLabel', data.ylabel_long(maskD) );
rotateticklabel(gca,45);
grid;
legend({ 'F1', 'F2', 'F3', 'F4'});
pause;

% idxD = find(data.C == 'c' | data.C == 'o'); % categorical and ordinal
% figure(2); subplot(3,1,1:2); hold off
% colors = {'k', 'b', 'r', 'c', 'm'};
% for dd=1:length(idxD)
%     d = idxD(dd);
%     ma = ~isnan(data.X(:,d));
%     U = unique(data.X(ma,d));
%     aux(size(Z,2)-1,length(U));
%     for f=2:size(Z,2) % for each feature        
%         x1 = data.X(Z(:,f) == 1,idxD);
%         x2 = data.X(Z(:,f) == 0,idxD);
%         aux(f-1,:) = hist(data.X(ma & (C == p),d),U);
%     end
%     aux = aux ./ repmat( sum(aux,2), 1, size(aux,2) );
% end

%% Plot per feature prob. of category activation
data.cat_labels{1} = data.cat_labels{1}';

idxs = find(data.C == 'c');

figure(3); hold off
leg_handlers = zeros(1,size(Z,2)); % number of lines to comment in legend
for k=1:size(Z,2) % for each feature
    ticklabels = [];
    xx = data.X(Z(:,k) == 1,:);
    idx_plot = 0;
    %subplot(5,1,1:4);
    for rr=1:length(idxs) % for each categorical dimension
        r = idxs(rr);
        V = xx(:,r);
        V(V == missing) = []; % remove missing for analysis
        h = hist(V,unique(data.X(~isnan(data.X(:,r)),r))) ./ length(V);
        p = plot(idx_plot+1:(idx_plot+length(h)), h,colors{k}, 'Linewidth', 2); hold on
        idx_plot = idx_plot+length(h);
        if isempty(data.cat_labels{r})
            lab = cell(length(h),1);
            for j=1:length(h)
                lab{j} = [data.ylabel_long{r}, ' ', num2str(j-1)];
            end
        else
            lab = data.cat_labels{r};
        end
        ticklabels = [ticklabels; lab];
        if (rr == 1)
            leg_handlers(k) = p;
        end
    end
end

legend(leg_handlers, {'F0', 'F1', 'F2', 'F3', 'F4'});
set( gca, 'XTick', 1:idx_plot)
set( gca,'XTickLabel', ticklabels );
rotateticklabel(gca,45);
grid;
pause;

%% Independent histograms for each identified pattern
for d=1:D % for each dimension
    ma = ~isnan(data.X(:,d));
    figure; hold off;
    if (data.C(d) == 'c') % plot histogram categories
        U = unique(data.X(ma,d));
        aux = zeros( size(patterns, 1),length(U));
        for p=1:size(patterns,1)
            aux(p,:) = hist(data.X(ma & (C == p),d),U);
        end
        aux = aux ./ repmat( sum(aux,2), 1, size(aux,2) );
        subplot(5,1,1:4);
        bar(unique(data.X(ma,d)), aux'); title(data.ylabel_long{d});
        set( gca,'XTickLabel', data.cat_labels{d}, 'FontSize', 14);
        rotateticklabel(gca,45);
        grid;
    else % plot box plots
        boxplot(data.X(ma,d),C(ma)); title(data.ylabel_long{d});
        grid;
    end
    pause;
end