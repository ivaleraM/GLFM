function [Xtrain,Xtest,Ytrain, Ytest] = splitData(X,Y,idx_Placebo,idx_Drug,perc,seed)
    % perc in training
    rng(seed)
    
    X_placebo = X(idx_Placebo,:);
    Y_placebo = Y(idx_Placebo,:);
    
    X_drug = X(idx_Drug,:);
    Y_drug = Y(idx_Drug,:);
    
    Ntrain_placebo = round(perc*size(X_placebo,1));
    Ntrain_drug = round(perc*size(X_drug,1));
    
    v_placebo = randperm(size(X_placebo,1));
    v_drug = randperm(size(X_drug,1));
    
    Xtrain = [ X_placebo(v_placebo(1:Ntrain_placebo),:); X_drug(v_drug(1:Ntrain_drug),:) ];
    Ytrain = [ Y_placebo(v_placebo(1:Ntrain_placebo),:); Y_drug(v_drug(1:Ntrain_drug),:) ];
    
    Xtest = [ X_placebo(v_placebo(1+Ntrain_placebo:end),:); X_drug(v_drug(1+Ntrain_drug:end),:) ];
    Ytest = [ Y_placebo(v_placebo(1+Ntrain_placebo:end),:); Y_drug(v_drug(1+Ntrain_drug:end),:) ];
    
    % finally, shuffle data
    v = randperm(size(Xtrain,1));
    Xtrain = Xtrain(v,:);
    
    v = randperm(size(Xtest,1));
    Xtest = Xtest(v,:);
    
    v = randperm(size(Ytrain,1));
    Ytrain = Ytrain(v,:);
    
    v = randperm(size(Ytest,1));
    Ytest = Ytest(v,:);
    
end