function plot_usa_map(data,Zn)

statesRef = data.cat_labels{1};
feature_level = zeros(length(statesRef),1);
for r=1:length(statesRef)
    idxs = find(data.X(:,1) == r); % list of counties corresponding to that state
    feature_level(r) = sum(Zn(idxs)) / length(idxs);
end

% statesRef: cell N*1
% feature_level: vector N*1

figure;
ax = usamap('all');
set(ax, 'Visible', 'off')
states = shaperead('usastatelo', 'UseGeoCoords', true);
for k = 1:numel(states)
    tmp = states(k).Name;
    idx = strncmp(tmp,statesRef,length(tmp));
    if isempty(idx)
        error('Unkown USA state');
    end
    states(k).Feature = feature_level(idx);
end
names = {states.Name};
indexHawaii = strcmp('Hawaii',names);
indexAlaska = strcmp('Alaska',names);
indexConus = 1:numel(states);
indexConus(indexHawaii|indexAlaska) = [];
surfaceColors = makesymbolspec('Polygon', {'Feature', [0 1], 'FaceColor', flip( autumn(numel(states)) ) });
%surfaceColors = makesymbolspec('Polygon', {'Feature', [min([states.Feature]) max([states.Feature])], 'FaceColor', flip( autumn(numel(states)) ) });

geoshow(ax(1), states(indexConus),  'SymbolSpec', surfaceColors);
geoshow(ax(2), states(indexAlaska),  'SymbolSpec', surfaceColors);
geoshow(ax(3), states(indexHawaii),  'SymbolSpec', surfaceColors);
colormap( flip( autumn(numel(states)) ) );
colorbar('Fontsize', 30);
%caxis( [min([states.Feature]) max([states.Feature])] )
caxis( [0 1] )
%stateColor = [0.5 1 0.5];
%geoshow(ax(1), states(indexConus),  'FaceColor', stateColor)
%geoshow(ax(2), states(indexAlaska), 'FaceColor', stateColor)
%geoshow(ax(3), states(indexHawaii), 'FaceColor', stateColor)

for k = 1:3
    setm(ax(k), 'Frame', 'off', 'Grid', 'off',...
        'ParallelLabel', 'off', 'MeridianLabel', 'off')
end