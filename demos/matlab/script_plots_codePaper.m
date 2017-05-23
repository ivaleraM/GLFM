%% script_generate_figures_codePaper
addpath(genpath('/home/melanie/Documents/UC3M/Workspace/additional_software/matlab2tikz-master'))
clear

%%
clear
%uiopen('/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/GLFM/demos/matlab/figs/transform2/yhist2.fig',1)
figure;
hold on;
load('./results/no_transform_counties.mat');
[h, xx] = hist(f_p_1(data.X(:,10), hidden.mu(10), hidden.w(10)), 20);
h = h ./ sum(h * (xx(2) - xx(1)));
bar(xx, h);
set(get(gca,'child'),'FaceColor',[1 1 1], 'FaceAlpha', 0, ...
    'EdgeColor',[0 102 255]/255);
hold on;

load('./results/counties_bias1_simId2_Niter10000_s2B1.00_alpha1.mat')
[h, xx] = hist(f_p_1(params.t_1{10}(data.X(:,10)), hidden.mu(10), hidden.w(10)), 20);
h = h ./ sum(h * (xx(2) - xx(1)));
bar(xx, h);
R = get(gca,'child');
set(R(1),'FaceColor',[1 1 1], 'FaceAlpha', 0, ...
    'EdgeColor',[0 128 0]/255);
hold on;

%%
load('./results/transform_counties.mat');
    
colors = [ 0 102 255; ...
            153 51  255; ...
            204 204 0; ...
            255 102  102; ...
            0   204 102];
    colors = colors ./ 255; %repmat(sum(colors,2),1,3);
    colors(3,:) = [0.9290    0.6940    0.1250];
    colors(5,:) = [0.4660    0.6740    0.1880];
    
    colors = colors([3 5 4 2 1],:);
    colors(4,:) = [0 255 255] ./255;
    
xx = -10:0.01:10;
V = [1 3];
dim = 10; % dimension to plot = white
hold on
for k=V % patterns to plot
    plot(xx,normpdf(xx,Zp(k,:)*hidden.B(dim,:,1)',hidden.s2Y(dim)), ...
        'Linewidth', 2, 'Color', colors(k,:));
    hold on
end
grid;
%legend({'T: 1  0  0  0', 'T: 0  1  0  0'})
saveas(gca,'./figs/transform2/yhist.fig')
xlim([-3 6])
legend({'h1', 'h2', '0  0  0  0', '1  0  0  0'})
legend('Location', 'Northwest')
%cleanfigure
%matlab2tikz('figs/transform2/yhist.tikz')

%%
load('./results/no_transform_counties.mat');
xx = -10:0.01:10;

% colors = [ 0 102 255; ... % azul
%             %153 51  255; ... % violeta
%             204 204 0; ... % verdamarillo
%             %255 102  102; ... % rojo
%             %0   204 102]; ...
%             ]/ 255;
        colors = [ 0 102 255; ...
            153 51  255; ...
            204 204 0; ...
            255 102  102; ...
            0   204 102];
    colors = colors ./ 255; %repmat(sum(colors,2),1,3);
    colors(3,:) = [0.9290    0.6940    0.1250];
    colors(5,:) = [0.4660    0.6740    0.1880];
    
    colors = colors([3 5 4 2 1],:);
    colors(4,:) = [0 255 255] ./255;
    
Zp = eye(4);
Zp = [ones(4,1),Zp];
    
V = [1 3]; %[1 3];
dim = 10; % dimension to plot = white
%figure(20);
for k=2:3 % patterns to plot
    plot(xx,normpdf(xx,Zp(k,:)*hidden.B(dim,:,1)',hidden.s2Y(dim)), ...
        'Linewidth', 2, 'Color', colors(V(k-1),:), 'LineStyle', '-.');
    hold on
end
xlim([-3 6])
%xlabel('y');
legend({'h1', 'h2', '0  0  0  0', '1  0  0  0'})

%%
% script_generate_plots_paper_counties;
% %load('./results/counties_bias1_simId2_Niter10000_s2B1.00_alpha1.mat');
% %[patterns, C, L] = get_feature_patterns_sorted(hidden.Z);
% figure(10)
% xlim([60 100])
% legend({'Empirical', '1  0  0  0', '0  1  0  0'})
% xlabel('x')
% saveas(gca,'./figs/transform2/xhist.fig')
% 
% cleanfigure
% matlab2tikz('figs/transform2/transform_white.tikz')

cleanfigure;
matlab2tikz('figs/transform2/yhist3.tikz')

%% NO TRANSFORM
uiopen('/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/GLFM/demos/matlab/figs/transform2/dim10_no_transform.fig',1)
R = get(gca,'child');
set(R(3),'FaceColor',[1 1 1], 'FaceAlpha', 0, ...
    'EdgeColor',[0 102 255]/255);
cleanfigure;
saveas(gca,'no_transform_white2.fig');
matlab2tikz('figs/transform2/no_transform_white2.tikz')

%% TRANSFORM
uiopen('/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/GLFM/demos/matlab/figs/transform2/transform_white.fig',1)
R = get(gca,'child');
set(R(3),'FaceColor',[1 1 1], 'FaceAlpha', 0, ...
    'EdgeColor',[0 128 0]/255);
saveas(gca,'transform_white2.fig');
matlab2tikz('figs/transform2/transform_white2.tikz')

