%% Script to create pcb database

clear
addpath('./aux/');
% load database and categorical columns
load_matrix_counties;
load_char_columns_counties;
load_header_counties;
header{8} = 'age6574'; header{9} = 'age75';
%header = import_header(file, 1, 1);

nameVars = {'county','state'};
idxs = [1,2];

% fill out categorical dimensions
info = cell(1,size(X,2));

for c=1:length(nameVars)
    str = sprintf('tmp = unique(%s);',nameVars{c});
    eval(str);
    for r=1:length(tmp)
        if strcmp(tmp{r},'') % empty row
            str = sprintf('X(strcmp(%s,tmp(r)),idxs(c)) = NaN;', nameVars{c});
        else
            str = sprintf('X(strcmp(%s,tmp(r)),idxs(c)) = r;', nameVars{c});
        end
        eval(str);
    end
    tmp(strcmp(tmp,'')) = [];
    info{idxs(c)} = tmp;
end

%bin_dims = [9,17]; % dimensions that are binary
%X(:,bin_dims) = X(:,bin_dims) + 1; % binary vars {0,1} should be {1,2}

% % remove undesidered dimensions (e.g., patient ids)
% data.xlabel = X(:,1);
% dims_to_remove = [1,18];
% X(:,dims_to_remove) = [];
% header(dims_to_remove) = [];
% info(dims_to_remove) = [];

data.X = X;
data.C = ['ccnnn','ngppn','pnppp','pppp'];
data.cat_labels = info;
data.ylabel = header;

%%
load('../../../datasets/mat/states.mat');

V = char(data.cat_labels{1});
for r=1:size(V,1)
    co = V(r,:);
    idx_sta = strncmp(co, states.code,length(co));
    co
    states.name{idx_sta}
    data.cat_labels{1}{r} = states.name{idx_sta};
end


save('/Users/melanie/Documents/GLFM/datasets/mat/counties.mat','data');