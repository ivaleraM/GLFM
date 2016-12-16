missing=-1;
p=0.1; %probability of missing
f_1=@(x,w) log(exp(w*x)-1);
df_1=@(x,w) w/(1-exp(-w*x));
logN=@(x) max(log(x),-30); %To avoid -Inf

%% Selecting missing data
randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));
load databases/Wine.mat
Nmiss=round(N*D*p);
miss=randperm(N*D);
miss=miss(1:Nmiss);
Xmiss=X;        % Observation matrix
Xmiss(miss)= missing; % Missing data are coded as missing
X(isnan(X))=missing;
[N, D]= size(X);
s2Y=1;   % Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
s2B=1;   % Variance of the Gaussian prior of the weigting matrices B
alpha=1; % Concentration parameter of the IBP
Nsim=5000; % Number of iterations for the gibbs sampler
bias = 0;
maxK= D;
%% Inference
Zini=double(rand(N,2)>0.8);
[Zest B Theta]= IBPsampler(Xmiss,C,Zini,bias,s2Y,s2B,alpha,Nsim,maxK,missing);


%% Compute test log-likelihood
XT=X;
ii=0;
for i=miss
    ii=ii+1;
    if (XT(i)~=-1 || ~isnan(XT(i)))
        d=ceil(i/N);
        n=mod(i,N);
        if (n==0)
            n=N;
        end
        j=j+1;
        Br=squeeze(B(d,:,1));
        if (C(d)=='g') 
            TLK(ii) = logN(normpdf(XT(i),Zest(:,n)'*Br',1));
        elseif (C(d)=='p' ) 
            TLK(ii) = logN(normpdf(f_1(XT(i),W(d)),Zest(:,n)'*Br',2))+logN(abs(df_1(XT(i),W(d))));
        elseif (C(d)=='c') 
           Br=squeeze(B(d,:,:));                    
           prob=zeros(1,R(d));
           for r=1:R(d)-1
               for r2=1:R(d)
                    if r2~=r
                        prob(r) =prob(r)+ logN(normcdf(Zest(:,n)'*(Br(:,r)-Br(:,r2)),0,1));
                    end
               end
           end 
           if 1-sum(exp(prob(1:R(d)-1)))<0
                prob(end)=-30;
                prob=logN(exp(prob)/sum(exp(prob)));
           else
                prob(end)=logN(1-sum(exp(prob(1:R(d)-1))));
           end
           TLK(ii) =prob(XT(i));

        elseif (C(d)=='o' ) 
            Br=squeeze(B(d,:,1));
            if XT(i)==1
                TLK(ii) = logN(normcdf(Theta(d,XT(i))-Zest(:,n)'*Br',0,1));
            elseif XT(i)==R(d)
                TLK(ii) = logN(1- normcdf(Theta(d,XT(i)-1)-Zest(:,n)'*Br',0,1));
            else
                TLK(ii) = logN(normcdf(Theta(d,XT(i))-Zest(:,n)'*Br',0,1)-normcdf(Theta(d,XT(i)-1)-Zest(:,n)'*Br',0,1));
            end
        elseif (C(d)=='n')     
             Br=squeeze(B(d,:,1));
            if XT(i)==0
                TLK(ii) = logN(normcdf(f_1(XT(i)+1,W(d))-Zest(:,n)'*Br',0,1));
            else
                TLK(ii) = logN(normcdf(f_1(XT(i)+1,W(d))-Zest(:,n)'*Br',0,1)-normcdf(f_1(XT(i),W(d))-Zest(:,n)'*Br',0,1));
            end
        end
    end
end  


