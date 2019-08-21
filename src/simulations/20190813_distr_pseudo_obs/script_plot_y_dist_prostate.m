%% script_generate_figures_codePaper
addpath(genpath('/home/melanie/Documents/UC3M/Workspace/additional_software/matlab2tikz-master'))
clear

%% PLOT PSEUDO-OBS FOR ALL OBSERVATIONS
savepath = './figs/resub3/prostate/';
%file_to_load = './results/prostateRed_bias1_simId2_Niter10000_s2Y0.00_s2B1.00_alpha1.00.mat';
file_to_load = '/home/melanie/results/glfm/pseudo-obs/prostateRed_bias1_simId1_Niter10000_s2B1.00_alpha1.00.mat';


load(file_to_load);
n_dims = size(data.transformed.X);
dimensions = [1,2,3,4,5]; %1:5;
hidden.mu = zeros(length(hidden.mu));
for kk=dimensions
    figure(kk);
    hold off;
    if isempty(params.t_1{kk})
        if (data.transformed.C(kk) == 'p')
            y_vec = f_p_1(data.transformed.X(:,kk), hidden.mu(kk), hidden.w(kk));
        elseif  (data.transformed.C(kk) == 'n')
            %[h, xx] = hist(data.transformed.X(:,kk), 200);
            %% scale empirical pdf
            %h = h .* df_p_1(xx,hidden.mu(kk),hidden.w(kk));
            y_vec = f_n_1(data.transformed.X(:,kk), hidden.mu(kk), hidden.w(kk));
        elseif (data.transformed.C(kk) == 'c')
            B_vec = squeeze(hidden.B(kk,:,:)); % K x R
            y_vec = zeros(size(data.transformed.X,1),1);
            for n=1:size(data.transformed.X,1)
                y_vec(n) = randn*sqrt(hidden.s2Y(kk)) + hidden.Z(n,:)* B_vec(:,data.transformed.X(n,kk));
            end
        elseif (data.transformed.C(kk) == 'o')
            B_vec = squeeze(hidden.B(kk,:,1))'; % K x R
            %y_vec = zeros(size(data.transformed.X,1),1);
            %for n=1:size(data.transformed.X,1)
            %    y_vec(n) = randn*sqrt(hidden.s2Y(kk)) + hidden.Z(n,:)* B_vec(:,data.transformed.X(n,kk));
            %end
            y_vec = randn(N,1)*sqrt(hidden.s2Y(kk)) + hidden.Z* B_vec;
        end
    else
        y_vec = f_p_1(data.transformed.X(:,kk), hidden.mu(kk), hidden.w(kk));
        %if (data.transformed.C(kk) == 'p')
        %    [h, xx] = hist(f_p_1(params.t_1{kk}(data.transformed.X(:,kk)), hidden.mu(kk), hidden.w(kk)), 20);
        %elseif  (data.transformed.C(kk) == 'n')
        %    [h, xx] = hist(f_n_1(params.t_1{kk}(data.transformed.X(:,kk)), hidden.mu(kk), hidden.w(kk)), 20);
        %end
    end
    %histogram(y_vec,'Normalization','probability');
    [h,xx] = hist(y_vec,20);
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
            plot(xx,nk(k)*normpdf(xx,patterns(k,:)*hidden.B(dim,:,1)',hidden.s2Y(dim)), ...
                'Linewidth', 2, 'Color', colors(k,:));
        end
        total_pdf = total_pdf + nk(k)*normpdf(xx,patterns(k,:)*hidden.B(dim,:,1)',hidden.s2Y(dim));
        hold on
    end
    grid;
    plot(xx,total_pdf, 'Linewidth', 2, 'Color', 'r','LineStyle', '-.');
    %legend({'T: 1  0  0  0', 'T: 0  1  0  0'})
    %saveas(gca,strcat(savepath,['yhist-k=',str(kk),'.fig']))
    xlim([-3 6])
    legend('off')
    %legend({'Empirical', '0 0 0 0', '0 1 0 0', '1 0 0 0', '0 0 1 0', '0 0 0 1', 'Predicted'})
    %legend('Location', 'Northwest')
    title(data.transformed.ylabel(kk))
    %cleanfigure
    %matlab2tikz('figs/transform2/yhist.tikz')
        cleanfigure;
    matlab2tikz([savepath,'yhist-k=',num2str(kk),'.tikz'])
    
end
