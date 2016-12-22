%% Script to create pcb database

clear
% load database and categorical columns
file = '/home/melanie/Documents/UC3M/Workspace/GIBP_Isabel/simulations/data/csv_xls/pbc.csv';
X = importfile(file, 2, 419);
load_char_columns_pbc;
header = import_header_semicolon(file, 1, 1);

% sex -> 5
% spiders -> 8
% hepatom -> 9
% ascites -> 10
% drug -> 16
% edema -> 18
nameVars = {'sex','spiders', 'hepatom', 'ascites', ...
    'drug', 'edema'};
idxs = [5,8,9,10,16,18];

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

bin_dims = 17;
X(:,bin_dims) = X(:,bin_dims) + 1; % binary vars {0,1} should be {1,2}

data.X = X;
data.C = ['ppopc','npccc','ppnnn','cccn'];
data.cat_labels = info;
data.dlabel = header;

save('./mat/pbc.mat','data');