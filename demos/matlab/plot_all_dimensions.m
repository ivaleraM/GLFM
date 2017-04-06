function plot_all_dimensions(data, hidden, params, Zp, leg, varargin)

if (length(varargin) == 1)
    idxD = varargin{1};
else
    idxD = 1:size(data.X,2);
end

for d=idxD
    figure;
    %  subplot(2,1,1);
    [xd, pdf]     = IBPsampler_PDF(data, Zp, hidden, params, d);
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        h = bar(pdf');
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