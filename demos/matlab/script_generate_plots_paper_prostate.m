%% script generate plots paper

% load('./results/prostateRed_bias1_simId2_Niter10000_s2Y0.00_s2B1.00_alpha1.mat');
% idx_featD = 4;

load('./results_APR2/prostateRed_bias1_simId2_Niter10000_s2Y0.00_s2B1.00_alpha1.mat');
idx_featD = 3;

sum(hidden.Z)
    feat_toRemove = find(sum(hidden.Z) < N*0.03);
    hidden = remove_dims(hidden, feat_toRemove);
    sum(hidden.Z)
    [patterns, C] = get_feature_patterns(hidden.Z);
    
    Kest = size(hidden.B,2);
    Zp = eye(Kest);
    %Zp(3,1) = 1;
    %Zp = [Zp; 0 1 1];
    Zp(:,1) = 1; % bias active
    Zp = Zp(1:min(5,Kest),:);
    leg = {' Empirical', ' F0',' F1', ' F2', ' F3', ' F4', ' F5'};
    colors = [ 0 102 255; ...
            153 51  255; ...            
            204 204 0; ...
            255 102  102; ...
            0   204 102];
    colors = colors ./ 255; %repmat(sum(colors,2),1,3);
    colors(3,:) = [0.9290    0.6940    0.1250];
    colors(5,:) = [0.4660    0.6740    0.1880];
    
%     colors = [  0    0.4470    0.7410; ...
%         0.8500    0.3250    0.0980; ...
%         0.9290    0.6940    0.1250; ...
%         0.4940    0.1840    0.5560; ...
%         0.4660    0.6740    0.1880];
    
for k=4:5
figure(k); hold off;
[h xx] = hist(data.X(:,k),100);
%xx(1) = 0.01;
%h = histogram(data.X(:,k),100);
h = h ./ sum(h * (xx(2) - xx(1)));
bar(h);
set(get(gca,'child'),'FaceColor',[0.8784 0.8784 0.8784], ...
    'EdgeColor',[0.7529 0.7529 0.7529]);
hold on;
end
set(gca,'xscale','log')

    for d=1:size(data.X,2)
            figure(d);
      %  subplot(2,1,1);
        [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
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
            h = stem(xd, pdf');
            for k=1:length(colors)
                h(k).Color
                h(k).Color = colors(k,:);
            end
        else
            if ~isempty(params.t{d})
                h = semilogx(xd,pdf','Linewidth', 2);
            else
                h = plot(xd,pdf','Linewidth', 2);
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
        legend(leg);
        grid;
    end