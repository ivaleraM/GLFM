import numpy as np
import random

def matrix_completion(Xmiss, C, s2Y=1, s2B=1, alpha=1, Niter=50, missing=-1):
    """
    Function to complete missing values of a certain numpy 2dim array

    Input parameters:
         Xmiss : numpy array which should be completed, should be normalized (TODO: check)
                 Size [NxD] where N is the number of observations and D is the
                 number of dimensions. Here missing data should be introduced
                 as the numeric value indicated in "missing", and categorical
                 and ordinal data should take values in {1,2, ..., R}.
         C     : char array [1xD] specifying the input data type of each column
                 (dimension) of the observation matrix X. Here 'g' indicates
                 real variable, 'p' positive real variable, 'n' count data,
                 'o' ordinal data and 'c' categorical data.

         s2Y   : variance of the Gaussian prior on the auxiliary variables
                 (pseudo-observations) Y (TODO: detail how to change it)
         s2B   : variance of the Gaussian prior on the elements of the
                 weighting matrices (latent features) B (TODO: Give intuition)
         alpha : mass parameter for the Indian Buffet Process
         Niter : number of internal iterations for the Gibbs sampler within
                 the C code before return
         missing : integer value that should be understood as missing value
    Output paramaters:
        Xcompl : same numpy array as Xmiss but whose missing values have been
                 inferred and completed by the algorithm.
    """

N = Xmiss.shape[0]
D = Xmiss.shape[1]
Xmiss[np.isnan(Xmiss)] = missing



## Inference
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
