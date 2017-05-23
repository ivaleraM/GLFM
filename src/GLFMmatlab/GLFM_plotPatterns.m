function GLFM_plotPatterns(data, hidden, params, patterns, varargin)

if ~isfield(data,'ylabel_long')
    data.ylabel_long = data.ylabel;
end

p = inputParser;
defaultColors = [];
defaultStyles = [];
defaultLeg = computeLeg(patterns, params);
defaultIdxD = 1:size(data.X,2);

checkIdxD = @(x) prod(ismember(x,1:size(data.X,2)));

addParameter(p,'idxD', defaultIdxD, checkIdxD);

% if specified, colors should be of length D or length(idxD)
checkColors = @(x) size(x,1) == size(patterns,1) || isempty(x); 
checkLeg = @(x) length(x) == size(patterns,1);

addRequired(p,'data',@isstruct);
addRequired(p,'hidden',@isstruct);
addRequired(p,'params',@isstruct);
addRequired(p,'patterns',@isnumeric);

addParameter(p,'colors', defaultColors,checkColors);
addParameter(p,'styles', defaultStyles,checkColors);

addParameter(p,'leg', defaultLeg,checkLeg);

parse(p, data, hidden, params, patterns, varargin{:});

%%
leg = [{'Empirical'}, p.Results.leg];
colors = p.Results.colors;
H_bar = cell(1,length(p.Results.idxD)); % figure handlers
H_cont = cell(1,length(p.Results.idxD));

for dd=1:length(p.Results.idxD)
    d = p.Results.idxD(dd);
    figure(d); hold off;
    % plot empirical if 'g' | 'p' | 'n'
    if ( data.C(d) == 'g' || data.C(d) == 'p' || data.C(d) == 'n')
        [h xx] = hist(data.X(:,d),100);
        %xx(1) = 0.01;
        %h = histogram(data.X(:,k),100);
        h = h ./ sum(h * (xx(2) - xx(1)));
        bar(xx, h);
        H_bar{dd} = h;
        set(get(gca,'child'),'FaceColor',[0.8784 0.8784 0.8784], ...
            'EdgeColor',[0.7529 0.7529 0.7529]);
        hold on;
    end
    %  subplot(2,1,1);
    [xd, pdf] = GLFM_computePDF(data, patterns, hidden, params, d);
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        mask = ~isnan(data.X(:,d)) & (data.X(:,d) ~= params.missing);
        tmp = hist(data.X(mask,d), unique(data.X(mask,d)));
        tmp = tmp / sum(tmp);
        h = bar([tmp' pdf']);
        h(1).FaceColor = [0.8784 0.8784 0.8784];
        H_bar{dd} = h;
    elseif (data.C(d) == 'n')
        h = plot(xd, pdf', 'Linewidth', 2);
        H_cont{dd} = h;
    else
        if ~isempty(params.t{d})
            h = semilogx(xd,pdf', 'Linewidth', 2);
        else
            h = plot(xd,pdf', 'Linewidth', 2);
        end
        H_cont{dd} = h;
        if isempty(colors)
            colors = zeros(size(h,1),3);
            for k=1:size(h,1)
                colors(k,:) = h(k).Color;
            end
        end
    end
    if ~isempty(colors)
        for k=1:size(h,1)
            if (data.C(d) == 'c') || (data.C(d) == 'o')
                ll = k+1;
                h(ll).FaceColor = colors(ll,:);
            else
                ll = k;
            end
            h(ll).Color = colors(ll,:);
        end
    end
    if ~isempty(p.Results.styles)
        for k=1:size(h,1)
            if isempty(p.Results.styles{k})
                continue;
            end
            if (data.C(d) == 'c') || (data.C(d) == 'o')
                ll = k+1;
            else
                ll = k;
            end
            h(ll).LineStyle = p.Results.styles{k};
        end
    end
    title(data.ylabel_long{d});
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        set(gca,'XTickLabel',data.cat_labels{d});
        set(gca,'XTickLabelRotation',45);
    end
    legend(leg);
    grid;
    %pause;
end

% if no color by default, set bar colors from histograms to line colors
for dd=1:length(p.Results.idxD)
    d = p.Results.idxD(dd);
    figure(d);
    if (data.C(d) == 'c') || (data.C(d) == 'o')
        for k=1:size(h,1)
            H_bar{dd}(k+1).FaceColor = colors(k,:);
        end
    end
end

function leg = computeLeg(Zp, params)
    if params.bias
        leg = num2str(Zp(:,2:end));
    else
        leg = num2str(Zp);
    end        
    leg = mat2cell(leg, ones(size(Zp,1),1), size(leg,2))';
end

end