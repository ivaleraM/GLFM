
function V = plot_cont_all_feats(data, hidden, params, varargin)
% last row of patterns corresponds to reference pattern

% Plot all features in same plot
%X is 2*D array with sctive vector + ref. vector

figure; hold on;
colors = {'b', 'r', 'm', 'c', 'g', 'y', 'b--', 'r--', 'm--', 'c--', 'g--', 'y--'};

Kest = size(hidden.Z,2);

idx_toRemove = [1 size(data.X,2)]; % CAREFUL: remove state and turnout
%idx_toRemove = [1 4 14 15 16]; % specific for counties DB

if isempty(varargin)

leg = cell(1,Kest-params.bias);
for k = params.bias+1:Kest
    Zp = zeros(2,Kest);
    Zp(:,1) = 1; % bias
    Zp(1,k) = 1; % feature active
    
    X_F = IBPsampler_MAP(data.C, Zp, hidden, params);
    
    labels = data.ylabel;
    
    X_F(:,idx_toRemove) = [];
    labels(idx_toRemove) = [];
    
    V = X_F(1,:) ./ X_F(2,:);
    V(V < 0) = 0;
    warning('Due to numerical errors, log of negative --> 0');
    V = log2( V );
    plot(V, colors{k-1}, 'Linewidth', 2);
    xticks(1:size(X_F,2));
    xticklabels(labels);
    xtickangle(90);
    leg{k-1} = sprintf('F%d',k-1);
end
legend(leg);
grid;

elseif (length(varargin) == 1)
    patterns = varargin{1};
    H = cell(1,size(patterns,1)); % container for plot handlers
    for k=1:(size(patterns,1))
        Zp = zeros(1,Kest);
        Zp(1,:) = patterns(k,:);
        %Zp(2,:) = patterns(end,:); % reference
        
        X_F = IBPsampler_MAP(data.C, Zp, hidden, params);
        
        labels = data.ylabel;
        
        X_F(:,idx_toRemove) = [];
        labels(idx_toRemove) = [];
    
        V = X_F(1,:) ./ mean_not_isnan(data.Xtrue(:,2:(end-1)),0); % X_F(2,:);
        V(V < 0) = 0;
        warning('Due to numerical errors, log of negative --> 0');
        V = log2( V );
        H{k} = plot(V, colors{k}, 'Linewidth', 2);
        xticks(1:size(X_F,2));
        xticklabels(labels);
        xtickangle(90);
        %leg{k} = sprintf('F%d',k-1);
    end
    legend({'(1 0 0)','(1 0 1)','(1 1 0)','(1 1 1)'});
    legend('Location','south');
    grid;
else
    error('Incorrect number of parameters');
end

% change color
colors2 = [ 0 102 255; ...
    153 51  255; ...
    204 204 0; ...
    255 102  102];% ...
%            0   204 102];
colors2 = colors2 ./ 255; %repmat(sum(colors,2),1,3);
colors2(3,:) = [0.9290    0.6940    0.1250];

for k=1:length(colors2)
    H{k}.Color = colors2(k,:);
end
ylim([-5 1]);
    
%plot([1 size(X_F,2)], [0 0], 'k', 'Linewidth', 2);