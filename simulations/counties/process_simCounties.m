%% Data Exploratory Anaylisis for Prostate DB
%addpath(genpath('./aux/'));
addpath(genpath('./'));

load('../results_old/tmp_counties.mat');
%load('../results/counties_Niter5000_1_s2Y1.00_s2u1.00_alpha1.00.mat');

fpos_1_handler = @(x,w) fpos_1(x,w);
dfpos_1_handler = @(x,w) dfpos_1(x,w);
fpos_handler = @(y, w) fpos(y,w);

data.ylabel_long = {'County', 'State', 'Msa', 'Pmsa', ...
    '1992 pop per 1990 miles^2', '1990 population', ...
    'perc population change 1980-1992', 'perc age 65-74, 1990', ...
    'perc age >= 75, 1990', 'serious crimes per 100,000 1991', ...
    'perc with bachelors degree or higher of those age>=25', ...
    'median family income, 1989 dollars', 'farm population, perc of total, 1990', ...
    'perc votes cast for democratic president', 'perc votes cast for republican president', ...
    'perc votes cast for Ross Perot', 'perc white, 1990', 'perc black, 1990', ...
    '1992 votes for president / 1990 pop x 100'};

Z = Zest';

th = 0.05; % remove features whose perc of observations is below th
idxs_to_delete = find( (sum(Z) > size(data.X,1)*(1-th)) | (sum(Z) < size(data.X,1)*th) );
Z(:,idxs_to_delete(2:end)) = [];
B(:,idxs_to_delete(2:end),:) = [];

%Z(:,2) = not(Z(:,2);

%[patterns, C] = get_feature_patterns([drug_identifier, Z]);
%print_patterns(1,data.X,patterns,C);
%Z(:,2) = [];
[patterns, C] = get_feature_patterns(Z);
idx = 14; % democratic vote
print_patterns(1,data.X(~isnan(data.X(:,idx)),idx),patterns,C,data.ylabel_long{idx});

%G1 = (C < 5);
%G2 = (C == 5) | (C == 6) | (C == 7);

%%

%% Plot theoretical distribution for each dimension
patterns = [ones(7,1), [zeros(1,6); eye(6)]];
relevantPatterns = 1:6;
leg = {'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7'};
s2u = 1;
colors = {'k--', 'g', 'b--', 'b', 'r--', 'r', 'm--', 'm'};
for d=5:size(data.X,2) % for each dimension
    figure(d); subplot(5,1,1:4); hold off;
    mask_noMissing = ~isnan( data.X(:,d) ) & data.X(:,d) ~= missing;
    if ((data.C(d) == 'n') || (data.C(d) == 'c') || (data.C(d) == 'o'))   && (min(data.X(:,d)) > 1)
        offset = min(data.X(mask_noMissing,d)) + 1;
    elseif (data.C(d) == 'p') && (min(data.X(:,d)) > 0)
        offset = min(data.X(mask_noMissing,d)) + 10^-6;
    end
    if (data.C(d) == 'c')
        R = length(data.cat_labels{d});
        pdfV = zeros(length(relevantPatterns),R);
    end
    for p=1:length(relevantPatterns) % for each feature activation pattern to compare
        pp = relevantPatterns(p);
        switch data.C(d)
            case 'g',
                xx = min(data.X(:,d)):0.01:max(data.X(:,d));
                Zn = patterns(pp,:);
                Bd = B(d,:,1)';
                pdf = pdf_real(xx, Zn,Bd,s2Y,s2u);
                plot(xx,pdf, colors{p}, 'linewidth', 2); hold on;
            case 'p',
                %xx = min(Xmiss(Xmiss(:,d)~=missing,d)):0.01:max(Xmiss(Xmiss(:,d)~=missing,d));
                xx = min(data.X(:,d)):0.1:max(data.X(:,d));
                Zn = patterns(pp,:);
                Bd = B(d,:,1)';
                w = 2 ./ max(data.X(:,d));
                pdf = pdf_pos(xx - offset,Zn,Bd,w,s2Y,s2u, fpos_1_handler, dfpos_1_handler);
                plot(xx,pdf, colors{p}, 'linewidth', 2); hold on;
            case 'c',
                Zn = patterns(pp,:);
                Bd = squeeze(B(d,:,:)); % TODO: Review that size = [K*R]
                pdf = pdf_cat(Zn,Bd,s2u,R);
                pdfV(p,:) = pdf;
            case 'o',
            case 'n',
                %xx = min(data.X(:,d)):1:max(data.X(:,d));
                xx = linspace(min(data.X(:,d)), max(data.X(:,d)),5000 );
                Zn = patterns(pp,:);
                Bd = B(d,:,1)';
                w = 2 ./ max(data.X(:,d));
                pdf = pdf_count(xx - offset,Zn,Bd,w,s2Y, fpos_1_handler);
                semilogx(xx,pdf, colors{p}, 'linewidth', 2); hold on;
            otherwise
                error('Unknown type of variable');
        end
    end
    if (data.C(d) == 'c')
        bar(1:R,pdfV', 'grouped');
%         set(gca,'XTick',[])
%         set(gca,'XTicklabel',[])
%         set( gca, 'XTick', 1:p)
        set( gca,'XTickLabel', data.cat_labels{d} );
        rotateticklabel(gca,45);
    end
    title(data.ylabel_long{d});
    legend(leg);
    grid;
end


%%

maskD = (data.C ~= 'c' & data.C ~= 'o'); % all dimensions excluding categorical and ordinal
figure(1); subplot(3,1,1:2); hold off
colors = {'k--', 'b', 'r', 'c', 'm', 'k', 'b--'};
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
legend({ 'F1', 'F2', 'F3', 'F4', 'F5', 'F6'}); % TODO: Generalize
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

% %% Plot per feature prob. of category activation
% data.cat_labels{1} = data.cat_labels{1}';
% 
% idxs = find(data.C == 'c');
% 
% figure(3); hold off
% leg_handlers = zeros(1,size(Z,2)); % number of lines to comment in legend
% for k=1:size(Z,2) % for each feature
%     ticklabels = [];
%     xx = data.X(Z(:,k) == 1,:);
%     idx_plot = 0;
%     %subplot(5,1,1:4);
%     for rr=1:length(idxs) % for each categorical dimension
%         r = idxs(rr);
%         V = xx(:,r);
%         V(V == missing) = []; % remove missing for analysis
%         h = hist(V,unique(data.X(~isnan(data.X(:,r)),r))) ./ length(V);
%         p = plot(idx_plot+1:(idx_plot+length(h)), h,colors{k}, 'Linewidth', 2); hold on
%         idx_plot = idx_plot+length(h);
%         if isempty(data.cat_labels{r})
%             lab = cell(length(h),1);
%             for j=1:length(h)
%                 lab{j} = [data.ylabel_long{r}, ' ', num2str(j-1)];
%             end
%         else
%             lab = data.cat_labels{r};
%         end
%         ticklabels = [ticklabels; lab];
%         if (rr == 1)
%             leg_handlers(k) = p;
%         end
%     end
% end
% 
% legend(leg_handlers, {'F0', 'F1', 'F2', 'F3', 'F4'});
% set( gca, 'XTick', 1:idx_plot)
% set( gca,'XTickLabel', ticklabels );
% rotateticklabel(gca,45);
% grid;
% pause;

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

%% Independent histograms for each identified pattern
for d=1:D % for each dimension
    ma = data.X(~isnan(data.X(:,d)),d);
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