BASEDIR_ALL = {'BioDeg','German','Internet','Nesarclighter','Wine'};
FILEN_ALL = {'bioDeg2','german','Internet2','Nesarc','Wine'};
DZ_ALL = {'10','20','50'};

%res_path = '/home/melanie/results/glfm/';
res_path = '/n/scratchlfs/doshi-velez_lab/melanie/results/glfm/';

for zz=1:length(DZ_ALL)
    Dz = DZ_ALL{zz};
    for bb=5:length(BASEDIR_ALL)
    clearvars -except BASEDIR_ALL bb FILEN_ALL DZ_ALL zz Dz res_path
    BASEDIR1 = BASEDIR_ALL{bb};
    FILEN1 = FILEN_ALL{bb};
    logN=@(x) max(log(x),-30);

    load([res_path, 'Dz=' Dz '/SimCluster' BASEDIR1 '/' (FILEN1) '_missP10_it1.mat']);
    try
        [N,D]=size(res.X_true);
    catch
        disp('I will pause now');
        who
        res
        pause;
    end
    Nit=20;
    Err=zeros(9,Nit);
    Err_D=zeros(D,9,Nit);
    nD=zeros(D,9,Nit);

    for Pmissing=[10:10:90]
        p=Pmissing/10;
        for it=1:Nit

            load([res_path, 'Dz=' Dz '/SimCluster' BASEDIR1 '/' (FILEN1)  '_missP' num2str(Pmissing) '_it' num2str(it) '.mat']);

            MU=res.X_pred;
            XT = res.X_true;
            miss=find(isnan(res.X_miss))';
            j=0;

            for i=miss
               if ~isnan(XT(i))
                d=ceil(i/N);
                j=j+1;
                nD(d,p,it)=nD(d,p,it)+1;
                if (res.C(d)=='g' )
                    Err_D(d,p,it)= Err_D(d,p,it)+ (XT(i) - MU(i))^2;
                elseif (res.C(d)=='p')
                    Err_D(d,p,it)= Err_D(d,p,it)+ (XT(i) - MU(i))^2;
                elseif (res.C(d)=='n' )
                    Err_D(d,p,it)= Err_D(d,p,it)+ (XT(i) - round(MU(i)))^2;
                elseif (res.C(d)=='c')
                    Err_D(d,p,it)= Err_D(d,p,it)+ (XT(i) ~= round(MU(i)));
                elseif (res.C(d)=='o')
                    Err_D(d,p,it)= Err_D(d,p,it)+ abs(XT(i) - round(MU(i)))/range(XT(:,d));
                end
               end
            end 
            Err(p,it)= 0;
            for d=1:D 
                if (res.C(d)=='g' || res.C(d)=='p' || res.C(d)=='n' ) 
                    Err_D(d,p,it) = sqrt(Err_D(d,p,it)/nD(d,p,it))/range(XT(:,d));
                    Err(p,it)=Err(p,it)+ Err_D(d,p,it);
                else
                    Err_D(d,p,it) = Err_D(d,p,it)/nD(d,p,it);
                    Err(p,it)=Err(p,it)+ Err_D(d,p,it);
                end
            end
            Err(p,it)=Err(p,it)/D;
            end
    end

    fileN= ['Dz=' Dz '_' 'Err_MDFA_' BASEDIR1 '.mat']; save(fileN, 'Err','Err_D')
    end
end
