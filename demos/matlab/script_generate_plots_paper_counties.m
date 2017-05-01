clear
close all
addpath(genpath('/home/melanie/Documents/UC3M/Workspace/additional_software/matlab2tikz-master'));

%load('./results/counties_bias0_simId1_Niter10000_s2B1.00_alpha1.mat');
%load('./results/counties_bias1_simId1_Niter10000_s2B1.00_alpha1.mat')
load('./results/counties_bias1_simId2_Niter10000_s2B1.00_alpha1.mat')
savepath = './figs/counties_bias1_mode2/';

%%
sum(hidden.Z)
feat_toRemove = find(sum(hidden.Z) < N*0.03);
hidden = remove_dims(hidden, feat_toRemove);
hidden = sort_hidden(hidden);
sum(hidden.Z)

idxF = 4; % 4; % 3;
%idxI = 4; % 5; % 4;

[patterns, C, L] = get_feature_patterns_sorted(hidden.Z);
[patterns, C, L] = get_feature_patterns_sorted(hidden.Z(:,1:idxF));
%%
idxD = 2:size(data.X,2);

Zp = eye(idxF);
Zp(:,1) = 1;
Zp = [Zp; 1 1 0 1]; %1 1 1 0];
%Zp = [Zp; 1 1 1];

Zp = patterns(L > 240,:);
leg = num2str(Zp(:,2:end));
leg = mat2cell(leg, ones(size(Zp,1),1), size(leg,2))';
Zp = [Zp zeros(size(Zp,1), size(hidden.Z,2) - idxF)];

colors = [ 0 102 255; ...
    153 51  255; ...
    204 204 0; ...
    255 102  102; ...
    0   204 102;
    255 51 255];
colors = colors ./ 255; %repmat(sum(colors,2),1,3);
colors(3,:) = [0.9290    0.6940    0.1250];
colors(5,:) = [0.4660    0.6740    0.1880];

% change order of colors
colors = colors([3 5 4 2 1],:);
colors(4,:) = [0 255 255] ./255;
%colors = colors([1 2 4 3 5 6],:);
%colors = colors([5 4 1 3 6 2],:);

styles = {'-','-', '-', '-','-','--'};
plot_all_dimensions(data, hidden, params, Zp, leg, colors, styles, idxD);

for d=idxD
    figure(d);
    if (d == 10)
        ylim([0 0.2]);
    elseif (d == 6)
        xlim([0 20]);
        ylim([0 0.25]);
    elseif (d == 10)
        xlim([80 100]);
        ylim([0 0.28]);
    elseif ((d == 7) || (d == 8) || (d == 9))
        xlim([5 85]);
    end
    if (d == 2)
        set(gca,'xscale','log');
        legend off
    elseif (d == 9) % Perot
        legend('Location','Northeast');
    elseif (d==10)
        legend('Location','Northwest');
    else
        legend off;
    end
    title('')
    cleanfigure;
    matlab2tikz([savepath, sprintf('dim%d.tex', d)] );
    saveas(gca,[savepath, sprintf('dim%d.fig', d)] );
    figurapdf(7,7)
    print([savepath, sprintf('dim%d.pdf', d)],['-d','pdf'],'-r150');
end
%%
for k=1:5%size(patterns,1)
    pat = patterns(k,:);
    Zn = (C == k);
    
    if (sum(Zn) < (size(data.X,1)*0.01))
        continue; % only plot patterns with more than X% of number of obs.
    end
    plot_usa_map(data,Zn);
    if (k ~=3)
        colorbar off;
    end
    %title(sprintf('Activation of pattern (%s)',num2str(pat)));
    cleanfigure;
    matlab2tikz([savepath, sprintf('clusters_map%d.tex', k)]);
    saveas(gca,[savepath, sprintf('clusters_map%d.fig', k) ] );
    figurapdf(7,7);
    print([savepath, sprintf('clusters_map%d.pdf', k) ],['-d','pdf'],'-r150');
end

% for k=2:4
%     Zn = (hidden.Z(:,k) == 1);
%     plot_usa_map(data,Zn);
%     if (k ~=4)
%         colorbar off;
%     end
%     matlab2tikz([savepath, sprintf('features_map%d.tex', k)]);
%     saveas(gca,[savepath, sprintf('features_map%d.fig', k) ] );
%     figurapdf(7,7);
%     print([savepath, sprintf('features_map%d.pdf', k) ],['-d','pdf'],'-r150');
% end