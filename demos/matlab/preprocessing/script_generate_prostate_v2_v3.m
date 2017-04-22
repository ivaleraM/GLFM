%% script to improve prostate.mat DB definition

input_file = '../../../datasets/mat/prostate.mat';
load(input_file);

%% ADAPT INPUT DATA --> put bias
data.cat_labels{1} = {'3';'4'};

% change to ordinal variables

% data.C(3) = 'o'; % months to follow up
% tmp = num2str( unique(data.X(~isnan(data.X(:,3)),3)) );
% data.cat_labels{3} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));

% Type of activity: confined in bed, in bed less than 1/2 daytime, ...
data.C(7) = 'o';
tmp = (data.X(:,7) == 2);
data.X( data.X(:,7) == 3,7) = 2;
data.X(tmp,7) = 3;
tmp = data.cat_labels{7}(2);
data.cat_labels{7}(2) = data.cat_labels{7}(3);
data.cat_labels{7}(3) = tmp;
%tmp = data.cat_labels_long{7}(2);
%data.cat_labels_long{7}(2) = data.cat_labels_long{7}(3);
%ata.cat_labels_long{7}(3) = tmp;

% Index Stage (state of the patient)
data.C(14) = 'o';
tmp = num2str( unique(data.X(~isnan(data.X(:,14)),14)) );
data.cat_labels{14} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));

% ---------------------
% drug_identifier = data.X(:,2) > 0.5;
% % remove drug levels
% data.X(:,2) = [];
% data.C(2) = [];
% data.cat_labels(2) = [];
% data.ylabel(2) = [];
% data.ylabel_long(2) = [];

tmp = num2str(unique(data.X(:,2)));
data.cat_labels{2} = mat2cell(tmp,ones(size(tmp,1),1),size(tmp,2));
data.C(2) = 'o';

save('../../../datasets/mat/prostate_v2.mat', 'data');

%% FURTHER MANIPULATIONS (simplification of some discrete variables

% simplify Status
data.cat_labels{4};
% 1 -->      1  alive
% [2,3] -->  2  vascular
% 6 -->      3  prostatic cancer
% [7,8] -->  4  lung-related dead
% [4,5,9,10] 5  others
V = data.X(:,4);
V(data.X(:,4) == 3) = 2;
V(data.X(:,4) == 6) = 3;
%V(data.X(:,4) == 7 | data.X(:,4) == 8) = 4;
%V(data.X(:,4) == 4 | data.X(:,4) == 5 | data.X(:,4) == 9 | data.X(:,4) == 10) = 5;
V(data.X(:,4) == 7 | data.X(:,4) == 8 | ...
    data.X(:,4) == 4 | data.X(:,4) == 5 | data.X(:,4) == 9 | data.X(:,4) == 10) = 4;
%data.cat_labels{4} = {'alive', 'vascular', 'prostatic ca', 'lung-related', 'others'};
data.cat_labels{4} = {'alive', 'vascular', 'prostatic ca', 'others'};
data.X(:,4) = V;

% for drug level, fusion 0 and 0.2
mask = logical(data.X(:,2) < 0.5);
data.X(mask,2) = 0;
data.cat_labels{2}(2,:) = [];

save('../../../datasets/mat/prostate_v3.mat', 'data');
