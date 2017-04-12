function plot_all_dimensions(data, hidden, params, Zp, leg, varargin)

if (length(varargin) == 1)
    idxD = varargin{1};
else
    idxD = 1:size(data.X,2);
end

for d=idxD
    figure(d); hold off;
    % plot empirical if 'g' | 'p' | 'n'
    if ( data.C(d) == 'g' || data.C(d) == 'p' || data.C(d) == 'n')
        [h xx] = hist(data.X(:,d),100);
        %xx(1) = 0.01;
        %h = histogram(data.X(:,k),100);
        h = h ./ sum(h * (xx(2) - xx(1)));
        bar(h);
        set(get(gca,'child'),'FaceColor',[0.8784 0.8784 0.8784], ...
            'EdgeColor',[0.7529 0.7529 0.7529]);
        hold on;
    end
    %  subplot(2,1,1);
    [xd, pdf] = IBPsampler_PDF(data, Zp, hidden, params, d);
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        mask = ~isnan(data.X(:,d)) & (data.X(:,d) ~= params.missing);
        tmp = hist(data.X(mask,d), unique(data.X(mask,d)));
        tmp = tmp / sum(tmp);
        h = bar([tmp' pdf']);
        h(1).FaceColor = [0.8784 0.8784 0.8784];
    elseif (data.C(d) == 'n')
        plot(xd, pdf');
    else
        if ~isempty(params.t{d})
            semilogx(xd,pdf', 'Linewidth', 2);
        else
            plot(xd,pdf', 'Linewidth', 2);
        end
    end
    title(data.ylabel_long{d});
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        set(gca,'XTickLabel',data.cat_labels{d});
        set(gca,'XTickLabelRotation',45);
    end
    legend(leg);
    pause;
end