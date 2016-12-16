function Xcompl= matrix_completion(Xmiss, C,s2Y, s2B, alpha, Niter, missing)
%% Inputs
% X:     Observation matrix of size $N x D$ where N is number of observations and
%        D is the dimension of the observations. Here missing data should
%        be introduced as the numeric value indicated in "missing", and categorical and ordinal data should take
%        values in {1,2, ..., R}.
% C:     Char array of length D, indicating the type of variable in each
%        column of the observation matrix X. Here 'g' indicates Real
%        variable, 'p' positive real variable, 'n' count data, 'o' ordinal
%        data and 'c' categorical data.
% missingVal: Value in which the missing data are encoded
% s2Y:   Variance of the Gaussian prior on the auxiliary variables (pseudoo-observations) Y
% s2B:   Variance of the Gaussian prior of the weigting matrices B
% alpha: Concentration parameter of the IBP
% Niter:  Number of iterations for the gibbs sampler
% missing:Missing value

[N D] = size(Xmiss);
W=ones(1,D);
for d=1:D
   if (C(d)=='g') 
        R(d) = 1;
    elseif (C(d)=='p' ) 
        R(d) = 1;
        W(d) = 1/max(Xmiss(:,d));
    elseif (C(d)=='c') 
       R(d) = max(Xmiss(:,d));
    elseif (C(d)=='o' ) 
        R(d) = max(Xmiss(:,d));
    elseif (C(d)=='n')     
        R(d) = 1;
        W(d) = 1/max(Xmiss(:,d));
   end 
end
maxR=max(R);
%% Inference
X(isnan(X))=missing;
Zini=double(rand(N,2)>0.8);
[Kest Zest B Y Theta]= IBPsampler(Xmiss,C,R,W,maxR,Zini,s2Y,s2B,alpha,Niter, missing);


%% Compute test log-likelihood
Xcompl=Xmiss;
miss=find(Xmiss==missing)';
f_1=@(x,w) log(exp(w*x)-1);
f=@(y,w) log(exp(y)+1)/w;
for ii=miss
    if Xmiss(ii)==missing
        d=ceil(ii/N);
        n=mod(ii,N);
        if (n==0)
            n=N;
        end
        Br=squeeze(B(d,:,1));
        if (C(d)=='g') 
            Xcompl(ii) = Zest(:,n)'*Br';
        elseif (C(d)=='p' ) 
            Xcompl(ii) = f(Zest(:,n)'*Br',W(d));
        elseif (C(d)=='c') 
           Br=squeeze(B(d,:,:));                    
           prob=zeros(1,R(d));
           Y=zeros(1,R(d));
           for r=1:R(d)
               Y(r)= Zest(:,n)'*Br(:,r);
           end 
           [val, Xcompl(ii)] = max(Y);
        elseif (C(d)=='o' ) 
            Br=squeeze(B(d,:,1));
            Y=Zest(:,n)'*Br';
            idx=find(Theta(d,1:R(d))>=Y);
            Xcompl(ii) = idx(1);
        elseif (C(d)=='n')     
            Br=squeeze(B(d,:,1));
            Xcompl(ii) = floor(f(Zest(:,n)'*Br',W(d)));
        end
    end
end  

end
