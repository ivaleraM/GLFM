%% Script to create pcb database

clear
% load database and categorical columns
file = '/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/simulations/data/csv_xls/prostate.csv';
load_matrix_prostate;
load_char_columns_prostate;
header = import_header(file, 1, 1);

nameVars = {'status','pf','ekg'};
idxs = [5,8,12];

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

bin_dims = [9,17]; % dimensions that are binary
X(:,bin_dims) = X(:,bin_dims) + 1; % binary vars {0,1} should be {1,2}

% remove undesidered dimensions (e.g., patient ids)
data.xlabel = X(:,1);
dims_to_remove = [1,18];
X(:,dims_to_remove) = [];
header(dims_to_remove) = [];
info(dims_to_remove) = [];

data.X = X;
data.C = ['cpncn','nccnn','cpnnp','c'];
data.cat_labels = info;
data.ylabel = header;

save('prostate.mat','data');