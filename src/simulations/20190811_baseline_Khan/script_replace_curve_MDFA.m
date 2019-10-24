%% script to replace MDFA curve
clear
close all

path_figs = '/home/melanie/git/GLFM/src/simulations/20190811_baseline_Khan/figures_error/';
new_path_figs = '/home/melanie/git/GLFM/src/simulations/20190811_baseline_Khan/new_figs_err/';
path_errs = '/home/melanie/git/GLFM/src/simulations/20190811_baseline_Khan/Err/';
mkdir(new_path_figs);

DB_FIG_NAMES = {'Wine','Internet','German','Nesarc','BioDeg'};
K_VEC = [10,20];

for dd=4:4%length(DB_FIG_NAMES)
    db = DB_FIG_NAMES{dd};
    
    % validation
    mean_val = zeros(1,length(K_VEC));
    for kk=1:length(K_VEC)
        K = K_VEC(kk);
        load([path_errs,'Dz=',num2str(K),'_Err_MDFA_',db,'.mat']);
        mean_val(kk) = mean(Err(:));
    end
    [~,best_kk] = min(mean_val);
    K = K_VEC(best_kk);
    fprintf('DB=%s, BEST K=%d\n',db,K);
    load([path_errs,'Dz=',num2str(K),'_Err_MDFA_',db,'.mat']);
    
   % open old figure
   uiopen([path_figs,'Err_',db,'2.fig']);
   hold on;
   figure(1);
   errorbar(1:1:9,mean(Err,2),2*std(Err'),'linewidth',3);
   legend({'GLFM','SIBP','BPFM','MDFA','MDFA'});
   savefig([new_path_figs,'Err_',db,'2.fig']);
   close(1);
 
   uiopen([path_figs,'Err_',db,'_D.fig']);
   hold on;
   figure(1);
   M = mean(Err_D,3);
   plot(1:size(M,1),M(:,5),'linewidth',3,'marker','x','linestyle',':');
   legend({'GLFM','SIBP','BPFM','MDFA','MDFA'});
   savefig([new_path_figs,'Err_',db,'_D.fig']);
   close(1);
end
    