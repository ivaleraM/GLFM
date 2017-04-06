%% script generate plots paper

% load('./results/prostateRed_bias1_simId2_Niter10000_s2Y0.00_s2B1.00_alpha1.mat');
% idx_featD = 4;

load('./results/counties_bias1_alpha1_simId10_Niter10000_s2B1.00_func1.mat');
data.ylabel_long = data.ylabel;

sum(hidden.Z)
feat_toRemove = find(sum(hidden.Z) < N*0.10);
hidden = remove_dims(hidden, feat_toRemove);
sum(hidden.Z)
[patterns, CC] = get_feature_patterns(hidden.Z);

for k=1:size(patterns,1)
    pat = patterns(k,:);
    Zn = (CC == k);%hidden.Z(:,idxF);
    plot_usa_map(data,Zn);
    title(sprintf('Activation of pattern (%s)',num2str(pat)));
end

plot_cont_all_feats(data, hidden, params);
for k=1:4
    figure(k);
    cleanfigure;
        matlab2tikz(sprintf('counties_map%d.tex', k));
        saveas(gca,sprintf('counties_map%d.fig', k) );
end
figure(5);
cleanfigure;
        matlab2tikz(sprintf('counties_featEffect.tex'));
        saveas(gca,sprintf('counties_featEffect.fig') );
   %% PLOT IND. DIMENSIONS
   close all;
   data.ylabel_long = data.ylabel;
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(6,Kest),:);
    leg = {'F0','F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
    
    Zp = patterns;
    leg = {'Empirical', '(1 0 0)', '(1 0 1)', '(1 1 0)', '(1 1 1)'};
    
        colors = [ 0 102 255; ...
            153 51  255; ...            
            204 204 0; ...
            255 102  102];% ...
%            0   204 102];
    colors = colors ./ 255; %repmat(sum(colors,2),1,3);
    colors(3,:) = [0.9290    0.6940    0.1250];
%    colors(5,:) = [0.4660    0.6740    0.1880];
    
idxD = [ 5 7 8 9 10 12 13];
piV = hist(CC,unique(CC));
piV = piV ./sum(piV);

    % Plot Empirical
    maxV = zeros(1,length(idxD));
    for kk=1:length(idxD)
        k = idxD(kk);
        if ~isempty(params.t{k})
            V = params.t_1{k}(data.X(:,k));
        else
            V = data.X(:,k);
        end
        figure(k); hold off;
        [h xx] = hist(V,100);
        h = h ./ sum(h * (xx(2) - xx(1)));
        maxV(kk) = max(h);
        bar(xx, h);
        set(get(gca,'child'),'FaceColor',[0.8784 0.8784 0.8784], ...
            'EdgeColor',[0.7529 0.7529 0.7529]);
        grid;
        hold on;
    end
    set(gca,'xscale','log')
maxV = ones(1,length(idxD));

    for kk=1:length(idxD)
        d = idxD(kk);
            figure(d);
        [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
        % normalize pdf to superpose empirical graphs
        %pdf = pdf ./ repmat( max(pdf,[],2), 1, size(pdf,2) );
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            mask = ~isnan(data.X(:,d));
            tmp = hist(data.X(mask,d), unique(data.X(mask,d)));
            tmp = tmp / sum(tmp);
            h = bar([tmp' pdf']);

            h(1).FaceColor = [0.8784 0.8784 0.8784];
            for k=1:length(colors)
                h(k+1).FaceColor = colors(k,:);
                if (k == idx_featD)
                    h(k+1).LineWidth = 2;
                    h(k+1).EdgeColor = [0.4016 0 0]; %[0.8 0 0]; %'red';
                end
            end
            
        elseif (data.C(d) == 'n')
            h = plot(xd, (pdf' .* repmat(piV, size(pdf,2),1) )  * maxV(kk),'Linewidth', 2);
            %h = plot(xd, pdf'  * maxV(kk),'Linewidth', 2);
            
            for k=1:length(colors)
                h(k).Color
                h(k).Color = colors(k,:);
            end
        else
            if ~isempty(params.t{d})
                %h = plot(xd,pdf' * maxV(kk),'Linewidth', 2);
                h = plot(xd,(pdf' .* repmat(piV, size(pdf,2),1) )  * maxV(kk),'Linewidth', 2);
            else
                %h = plot(xd,pdf'  * maxV(kk),'Linewidth', 2);
                h = plot(xd,(pdf' .* repmat(piV, size(pdf,2),1) )  * maxV(kk),'Linewidth', 2);
            end
            for k=1:length(colors)
                h(k).Color
                h(k).Color = colors(k,:);
            end
        end

        title(data.ylabel_long{d});
        if (data.C(d) == 'c') || (data.C(d) == 'o')
            set(gca,'XTickLabel',data.cat_labels{d});
            set(gca,'XTickLabelRotation',45);
        end
        grid;
        legend(leg);
        grid;
        cleanfigure;
        matlab2tikz(sprintf('counties_dim%d.tex', d));
        saveas(gca,sprintf('counties_dim%d.fig', d) );
    end