% function to run Khan baseline on different datasets
% Pmissing: percentage of missings
% mask_seed: random seed for missing mask (note: name is misleading!)
% Nsim: number of iterations to run
function run_baseline_mixed_FA(Pmissing, mask_seed, inCluster, db_id)
    tic;
    % load dataset and format accordingly
    if (db_id < 0 || db_id > 5)
        error('Incorrect db id')
    end
    if (inCluster == 0)
        BASEDIR = '/home/melanie/datasets/glfm/';
        RESDIR = '/home/melanie/results/glfm/';
        GLFMDIR = '/home/melanie/git/GLFM/';
    else
        BASEDIR = '/n/home11/mfernandezpradier/datasets/glfm/';
        RESDIR = '/n/home11/mfernandezpradier/results/glfm/';
        GLFMDIR = '/n/home11/mfernandezpradier/git/GLFM/';
    end
    addpath(genpath(pwd));
    addpath(genpath(fullfile(GLFMDIR,'/src/baselines/2010-NIPS-FA-code/')));
    
    db_paths = {'SimClusterBioDeg','SimClusterGerman','SimClusterInternet',...
        'SimClusterNesarclighter', 'SimClusterWine'};
    db_files = {'bioDeg2','german','Internet2','Nesarc','Wine'};
    datapath = fullfile(BASEDIR, db_paths{db_id});
    dbpath =  fullfile(BASEDIR, db_paths{db_id}, [db_files{db_id},'.mat']);
    fn = sprintf('missP%d_it%d.mat', Pmissing, mask_seed);
    misspath = fullfile(BASEDIR, db_paths{db_id}, fn);
    load(dbpath)
    load(misspath)
    % correct input format
    data.X_true = X;
    data.C = C;
    data.miss = miss;
    data.idxs_continuous = (data.C == 'p' | C == 'g');
    data.idxs_discrete = not(data.C == 'p' | C == 'g');
    data.continuous = data.X_true(:,data.idxs_continuous)';
    data.discrete = data.X_true(:,data.idxs_discrete)';
    Dcat = size(data.discrete,1);
    data.nClass = zeros(Dcat,1);
    for m=1:Dcat
        data.nClass(m) = length(unique(data.discrete(m,:)));
    end
    
    % standardize data (Important! Before missing mask!)
    data = standardize(data);    
    
    % apply missing mask
    Xmiss = data.X;
    Xmiss(data.miss) = nan;
    data.X_miss = Xmiss;
    data.continuous = Xmiss(:,data.idxs_continuous)';
    data.discrete = Xmiss(:,data.idxs_discrete)';
    
    % run mixed_FA_miss
    Dz = 5;
    res = wrapper_mixed_FA(data,Dz);
    
    % revert data standardization
    res.X_pred_mixedFA_standardized = standardize_back(res.X_pred_mixedFA,data);
    
    % change output format
    X_pred_ND = zeros(size(data.X_true));
    X_pred_ND(:,data.idxs_continuous) = res.X_pred_mixedFA_standardized.continuous';
    X_pred_ND(:,data.idxs_discrete) = res.X_pred_mixedFA_standardized.discrete';
    R = zeros(1,D);
    R(data.idxs_discrete) = data.nClass;
    
    % save results
    res.X_pred = X_pred_ND;
    res.X_true = data.X_true;
    res.X_miss = data.X_miss;
    res.C = data.C;
    res.R = R;
    toc
    res.err = compute_glfm_errors(res.X_true,res.X_pred,res.X_miss,res.C, res.R);
    respath = fullfile(RESDIR,db_paths{db_id});
    mkdir(respath)
    save(fullfile(respath,[db_files{db_id} '_' fn]),'res');

end