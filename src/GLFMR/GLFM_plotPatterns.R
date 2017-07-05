GLFM_plotPatterns<-function(data,hidden,params,Zp,varargin){
  
D<-dim(hidden$B[[1]])[1]
#print(D)
#for(dd in 1:length(p$Results$idxD)){
  for(d in 1:D){
 # d<-p$Results$idxD[dd]
  if(data$C[d]=='g'|| data$C[d] == 'p' || data$C[d] == 'n'){
    counts <- table(data$X[,d])
    barplot(counts)
    print("Press return to continue")
  }
    #duda aqui
    pdf_val<-GLFM_computePDF(data,Zp,hidden,params,d)
    if(data$C[d]=='c'||data$C[d] == 'o' ){
      idxs_missing<-union(which(data$X[,d] == params$missing),which(is.nan(data$X[,d])))
      idxs_not_missing<-setdiff(1:(dim(data$X)[2]),idxs_missing)
      #duda aqui
    }
    else if(data$C[d] == 'n'){
      h <-plot(pdf_val$xd)
      
    }
    else{
      h <-plot(pdf_val$xd)
    }
  }
}
  #To be completed


