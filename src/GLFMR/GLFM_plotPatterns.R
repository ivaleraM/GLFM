GLFM_plotPatterns<-function(data,hidden,params,varargin){
  
for(dd in 1:length(p$Results$idxD)){
  d<-p$Results$idxD[dd]
  if(data$C[d]=='g'||data$C[d] == 'p' || data$C[d] == 'n'){
    counts <- table(data$X[,dd])
    barplot(counts)
  }
    #duda aqui
    pdf_val<-GLFM_computePDF(data,patterns,hidden,params,dd)
    if(data$C[d]=='c'||data$C[d] == 'o' ){
      idxs_missing<-union(which(data$X[,d] == params$missing),which(is.nan(data$X[,d])))
      idxs_not_missing<-setdiff(1:(dim(data$X)[2]),idx_missing)
      #duda aqui
    }
    else if(data$C[d] == 'n'){
      h <-plot(pdf_val$xd)
      
    }
    else{
      h <-plot(pdf_val$xd)
    }
}
  #To be completed
}

