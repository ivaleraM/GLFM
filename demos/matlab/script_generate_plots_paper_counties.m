clear
close all
load('./results/counties_bias0_simId1_Niter10000_s2Y0.00_s2B1.00_alpha1.mat');
addpath(genpath('/home/melanie/Documents/UC3M/Workspace/additional_software/matlab2tikz-master'));

sum(hidden.Z)
feat_toRemove = find(sum(hidden.Z) < N*0.03);
hidden = remove_dims(hidden, feat_toRemove);
hidden = sort_hidden(hidden);
sum(hidden.Z)

idxF = 3; % 4; % 3;
idxI = 4; % 5; % 4;

[patterns, C] = get_feature_patterns(hidden.Z(:,1:idxF));

idxD = 2:size(data.X,2);

Zp = patterns(idxI:end,:);


%Zp = eye(size(hidden.Z,2));
leg = num2str(Zp);
leg = mat2cell(leg, ones(size(Zp,1),1), size(leg,2))';
Zp = [Zp zeros(size(Zp,1), size(hidden.Z,2) - idxF)];

plot_all_dimensions(data, hidden, params, Zp, leg, idxD);
for d=idxD
    figure(d);
    if (d == 2)
        set(gca,'xscale','log');
    elseif (d == 6)
        xlim([0 20]);
        ylim([0 0.5]);
    elseif (d == 10)
        xlim([80 100]);
        ylim([0 0.28]);
    elseif ((d == 7) || (d == 8) || (d == 9))
        xlim([5 85]);
    end
    cleanfigure;
    matlab2tikz(sprintf('./figs/counties_new/dim%d.tex', d));
    saveas(gca,sprintf('./figs/counties_new/dim%d.fig', d) );
    figurapdf(7,7)
    print(sprintf('./figs/counties_new/dim%d.pdf', d),['-d','pdf'],'-r300');
end

for k=1:size(patterns,1)
    pat = patterns(k,:);
    Zn = (C == k);
    %Zn = hidden.Z(:,k);
    
    if (sum(Zn) < (size(data.X,1)*0.01))
        continue; % only plot patterns with more than 5% of number of obs.
    end
    plot_usa_map(data,Zn);
    title(sprintf('Activation of pattern (%s)',num2str(pat)));
    cleanfigure;
    matlab2tikz(sprintf('./figs/counties_new/counties_map%d.tex', k));
    saveas(gca,sprintf('./figs/counties_new/counties_map%d.fig', k) );
    figurapdf(7,7)
    print(sprintf('./figs/counties_new/counties_map%d.pdf', k),['-d','pdf'],'-r300');
end
