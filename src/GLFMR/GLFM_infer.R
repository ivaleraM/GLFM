#' @param data is a list with X and C
#' @param varargin is a list of lists containing hidden variables and parameters 
#' @return A list with posterior draws

GLFM_infer<-function(data,varargin){
  source("init_default_params.R")
  varargin_size<-length(varargin)
  if(varargin_size==0){
    Z<-c()
    params <- init_default_params(data, c())
    print("Case 0")
  }
  else if(varargin_size<3)
  {
    Z<-varargin[1]
    readline("press return to continue")
    
    switch(varargin_size,params <- init_default_params(data, c()) ,params <- init_default_params(data, varargin[2]))
  }
  else
  {
    stop("Incorrect number of input parameters: should be 0, 1 or 2")
  }
  D <- dim(data$X)[2]
  N <- dim(data$X)[1]
  if(length(Z)==0)
  {
    m0<-matrix(0,N,2)
    Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
  }
  if(params$bias == 1 && length(params$bias)>0){
    Z <-cbind(rep(1,N),Z)
  }
  else if(params$bias>1  && length(params$bias)>0){
    stop("There is more than 1 bias specified, but the structure Z has not been initialized") 
  }
  # replace missing values
  idx_missing<-which(is.nan(data$X))
  data$X[idx_missing] <- params$missing
  # Change labels of categorical or ordinal data
  # Cambiamos las etiquetas de los datos categoricos o ordinales que no son datos
  # faltantes para que su respectiva categoria empiece en uno. y por eso se toma
  # el minimo, etc.
  idx_catdata<-which(data$C=='c')
  idx_orddata<-which(data$C=='o')
  set_of_both<-union(idx_catdata,idx_orddata)
  idx_not_missing<-setdiff(set_of_both,which(set_of_both%in%idx_missing))
  if(length(idx_not_missing)>0){
  V_offsets <- min( data$X[idx_not_missing] )
  data$X[idx_not_missing]<-data$X[idx_not_missing]-V_offsets+1
  }
  # eventually, apply external transform
  #if(t%in%params){
    #work in logarithm space
    #data.X(:,r) = params.t_1{r}(data.X(:,r)); # work in logarithm space better
    #data.C(r) = params.ext_dataType{r};
    # ---To be completed ---
  #}
 
  func_bit<-rep(1,dim(data$X)[2])
  
  # call .Rcpp wrapper function
  setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/")
  library(RcppGSLExample)
 hidden<-IBPsampler(t(data$X),unlist(data$C),t(Z),unlist(params$bias),func_bit,unlist(params$s2u),unlist(params$s2B),unlist(params$alpha),unlist(params$Niter),unlist(params$maxK), unlist(params$missing))
 #print(list(data$X, data$C, Z))
  print(list(t(data$X),Z))
# readline("press return to continue")
  #hidden<-IBPsampler(t(data$X),data$C,t(Z),params$bias,func_bit,params$s2u,params$s2B,
  #                  params$alpha,params$Niter,params$maxK,params$missing)
 R<-rep(1,D)
 hidden<-append(params,list("R"=R))
 setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
  return(hidden)
}
 
  
  
  
  