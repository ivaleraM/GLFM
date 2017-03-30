
function V = plot_cont_all_feats(data, hidden, params)
% Plot all features in same plot
%X is 2*D array with sctive vector + ref. vector

figure; hold on;
colors = {'b', 'r', 'm', 'c', 'g', 'y', 'b--', 'r--', 'm--', 'c--', 'g--', 'y--'};

Kest = size(hidden.Z,2);

idx_toRemove = [1 14]; % CAREFUL
%idx_toRemove = [1 4 14 15 16]; % specific for counties DB
leg = cell(1,Kest-params.bias);
for k = params.bias+1:Kest
    Zp = zeros(2,Kest);
    Zp(:,1) = 1; % bias
    Zp(1,k) = 1; % feature active
    
    X_F = IBPsampler_MAP(data.C, Zp, hidden);
    labels = data.ylabel;
     X_F(:,idx_toRemove) = [];
     labels(idx_toRemove) = [];
    
    V = log2( X_F(1,:) ./ X_F(2,:) );
    plot(V, colors{k-1}, 'Linewidth', 2);
    xticks(1:size(X_F,2));
    xticklabels(labels);
    xtickangle(45);
    grid;
    leg{k-1} = sprintf('F%d',k-1);
end
legend(leg);
grid;
%plot([1 size(X_F,2)], [0 0], 'k', 'Linewidth', 2);