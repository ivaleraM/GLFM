%% script_generate_figures_codePaper
addpath(genpath('/home/melanie/Documents/UC3M/Workspace/additional_software/matlab2tikz-master'))
clear

%% PLOT PSEUDO-OBS FOR ALL OBSERVATIONS
savepath = './figs/resub3/counties/';
%%%file_to_load = './results/no_transform_counties.mat';
%%file_to_load = './results/counties_bias1_simId2_Niter10000_s2B1.00_alpha1.mat';
%file_to_load = './results/transform_counties.mat';
file_to_load = '/home/melanie/results/glfm/pseudo-obs/counties_bias1_simId1_Niter10000_s2B1.00_alpha1.mat';
load(file_to_load);

%% normalize data
% Xmiss=data.X;        % Observation matrix
% 
% [N, D]= size(data.X);
% %s2Y=.5;    % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
% %s2u=.5;
% s2B=1;      % Variance of the Gaussian prior of the weigting matrices B
% %alpha=1;    % Concentration parameter of the IBP
% %Nsim=100; % Number of iterations for the gibbs sampler
% maxK= D;
% 
% Xmiss(isnan(Xmiss)) = missing;
% for d=1:D
%     % normalize data (such that it occupies interval (0,max(X(d))]
%     if ((data.C(d) == 'n') || (data.C(d) == 'c') || (data.C(d) == 'o'))   && (min(data.X(:,d)) > 1)
%         offset = min(Xmiss(Xmiss(:,d) ~= missing,d));
%         Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 1;
%     end
%     if (data.C(d) == 'p') && (min(data.X(:,d)) > 0)
%         offset = min(data.X(:,d));
%         Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) - offset + 10^-6;
%     end
% 
%     if ((data.C(d) == 'n' || (data.C(d) == 'c') || (data.C(d) == 'o')) && (min(data.X(:,d)) == 0))
%         Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 1;
%     elseif (data.C(d) == 'p') && (min(data.X(:,d)) == 0)
%         Xmiss(Xmiss(:,d) ~= missing,d) = Xmiss(Xmiss(:,d) ~= missing,d) + 10^-6;
%     end
%     
% end
%
%data.X = Xmiss;

%% Compute histogram of pseudo-observations
n_dims = size(data.transformed.X);
dimensions = 2:11;

for kk=dimensions
    %hidden.mu = min(data.transformed.X,[],1)';

    figure(kk);
    hold off;
    if isempty(params.t_1{kk})
        if (data.transformed.C(kk) == 'p')
            %f1x = f_p_1(data.transformed.X(:,kk), -3.6, hidden.w(kk));
            %f1x = f_p_1(data.transformed.X(:,kk), hidden.mu(kk), hidden.w(kk));
            f1x = f_p_1(data.X(:,kk), hidden.mu(kk), hidden.w(kk));
            s2_y_post = 1.0/(1.0/hidden.s2Y(kk) + 1.0/params.s2u);
            mu_y_post = ((hidden.Z *hidden.B(kk,:,1)')/hidden.s2Y(kk) + ...
               f1x/params.s2u) * s2_y_post;
            y_vec = randn(N,1)*sqrt(s2_y_post) + mu_y_post;
            y_vec = f1x;
        elseif  (data.transformed.C(kk) == 'n')
            %y_vec = f_n_1(data.transformed.X(:,kk), -3.6, hidden.w(kk));
            y_vec = f_n_1(data.transformed.X(:,kk), hidden.mu(kk,:), hidden.w(kk));
        elseif (data.transformed.C(kk) == 'c')
            continue;
        %    B_vec = squeeze(hidden.B(kk,:,:)); % K x R
        %    y_vec = zeros(size(data.transformed.X,1),1);
        %    for n=1:size(data.transformed.X,1)
        %        y_vec(n) = randn*sqrt(hidden.s2Y(kk)) + hidden.Z(n,:)* B_vec(:,data.transformed.X(n,kk));
        %    end
        %    [h, xx] = hist(y_vec,20);
        end
    else
        %[h, xx] = hist(f_p_1(params.t_1{kk}(data.transformed.X(:,kk)), hidden.mu(kk), hidden.w(kk)), 20);
        if (data.transformed.C(kk) == 'p')
            %y_vec = f_p_1(params.t_1{kk}(data.transformed.X(:,kk)),-0.7, hidden.w(kk));
            y_vec = f_p_1(data.transformed.X(:,kk),hidden.mu(kk), hidden.w(kk));
        elseif  (data.transformed.C(kk) == 'n')
            %y_vec = f_n_1(params.t_1{kk}(data.transformed.X(:,kk)),9.07, hidden.w(kk));
            y_vec = f_n_1(params.t_1{kk}(data.transformed.X(:,kk)),hidden.mu(kk,:), hidden.w(kk));
        end
    end
    %histogram(y_vec,10,'Normalization','probability');
    [h, xx] = hist(y_vec, 20);
    h = h ./ sum(h * (xx(2) - xx(1)));
    bar(xx, h);
    R = get(gca,'child');
    set(R(1),'FaceColor',[1 1 1], 'FaceAlpha', 0, ...
        'EdgeColor',[0 102 255]/255); %[0 128 0]/255);
    hold on;
    
    %%    
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
    [patterns, C, L] = get_feature_patterns_sorted(hidden.Z);
    nk = L / size(hidden.Z,1);
    V = [1 2 3 4 5];
    dim = kk; % dimension to plot = white
    hold on
    total_pdf = zeros(size(xx,1));
    for k=1:size(patterns,1) % patterns to plot
        if ismember(k,V)
            plot(xx,nk(k)*normpdf(xx,patterns(k,:)*hidden.B(dim,:,1)',sqrt(hidden.s2Y(dim))), ...
                'Linewidth', 2, 'Color', colors(k,:));
        end
        total_pdf = total_pdf + nk(k)*normpdf(xx,patterns(k,:)*hidden.B(dim,:,1)',sqrt(hidden.s2Y(dim)));
        hold on
    end
    
    grid;
    plot(xx,total_pdf, 'Linewidth', 2, 'Color', 'r','LineStyle', '-.');
    %legend({'T: 1  0  0  0', 'T: 0  1  0  0'})
    %saveas(gca,strcat(savepath,['yhist-k=',str(kk),'.fig']))
    xlim([-3 8])
    legend({'Empirical', '0  0  0', '1  0  0', '1  0  1', '0  1  0', '1  1  0', 'Predicted'})
    legend('Location', 'Northwest')
    %legend('off');
    title(data.transformed.ylabel(kk))
    %cleanfigure
    %matlab2tikz('figs/transform2/yhist.tikz')
    
    cleanfigure;
    matlab2tikz([savepath,'yhist-k=',num2str(kk),'.tikz'])
    
%     total_pdf2 = zeros(size(xx,1));
%     for n=1:200%size(data.transformed.X,1)
%         total_pdf2 = 1.0/size(xx,1)*normpdf(xx,mu_y_post(n),sqrt(s2_y_post));
%         plot(xx,total_pdf2, 'Linewidth', 2, 'Color', 'k','LineStyle', '-.');
%     end
    
end
