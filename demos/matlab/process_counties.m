%% process_counties

out_folder = './results/'; 
A = dir([out_folder, '*_Niter10000_*_func1.mat']);

for i=1:length(A)
    load([out_folder, A(i).name]);
    sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*0.1);
    hidden = remove_dims(hidden, feat_toRemove);
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    
    
    for k=1:size(patterns,1)
        pat = patterns(k,:);
        Zn = (C == k);%hidden.Z(:,idxF);
        plot_usa_map(data,Zn);
        title(sprintf('Activation of pattern (%s)',num2str(pat)));
    end
    
    
    for k=1:size(hidden.Z,2)
        Zn = hidden.Z(:,k);
        plot_usa_map(data,Zn);
        title(sprintf('Activation of feature %d',k));
    end
    
    plot_cont_all_feats(data, hidden, params);

   %% PLOT IND. DIMENSIONS
   
   data.ylabel_long = data.ylabel;
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(6,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
    
    Zp = patterns;
    leg = {'(1 0 0)', '(1 0 1)', '(1 1 0)', '(1 1 1)'};
    
    figure;
    for d=4:D
        subplot(2,1,1);
        [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            bar(pdf');
        elseif (data.C(d) == 'n')
            stem(xd, pdf');
        else
            plot(xd,pdf', 'Linewidth', 2);
        end
        title(data.ylabel_long{d});
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            set(gca,'XTickLabel',data.cat_labels{d});
            set(gca,'XTickLabelRotation',45);
        end
        legend(leg);
        subplot(2,1,2);
        if ~isempty(params.t{d})
            data.X(:,d) = params.t_1{d}(data.X(:,d));
        end
        hist(data.X(:,d),100); title('Empirical');
        pause;
    end
    
end
