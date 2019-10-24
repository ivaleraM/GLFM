clear
BASEDIR1 = ['Wine'];
logN=@(x) max(log(x),-30);

load(['SimCluster' BASEDIR1 '/' (BASEDIR1) '_missP10_it1.mat']);
[N,D]=size(res.X_true);
Nit=10;
Err=zeros(9,Nit);
Err_D=zeros(D,9,Nit);
nD=zeros(D,9,Nit);

for Pmissing=[10:10:90]
    p=Pmissing/10;  
    for it=1:Nit
        %try
        %load(['SimCluster' BASEDIR1 '/' lower(BASEDIR1) '.mat']);
                
        load(['SimCluster' BASEDIR1 '/' (BASEDIR1)  '_missP' num2str(Pmissing) '_it' num2str(it) '.mat']);%BMF
        %load(['IBP_General/' BASEDIR1 '/Results' results '/' BASEDIR1 'P' num2str(Pmissing) '_it' num2str(it) '.mat']);%'P'BASEDIR1
       
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
                %Err(p,it) =Err(p,it)+ abs(XT(i) - MU(i))/range(XT(:,d));
            end
           end
        end 
        Err(p,it)= 0;
%         dvec=ceil(miss/N);
%         [SND,aux] = histc(dvec,1:D);
%         SND(find(SND==0))=1;
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
%         catch err
%             disp('falta');
%     end
        
        end
        
end

fileN= ['Err_MDFA_' BASEDIR1 '.mat']; save(fileN, 'Err','Err_D')

    