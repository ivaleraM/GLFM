#var1<-c("HOLA")
# sapply(var1, tolower)
#varargin is a list
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
    switch(varargin_size,params <- init_default_params(data, c()) ,params <- init_default_params(data, unlist(varargin[2])))
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
  if(params$bias == 1){
    Z <-cbind(rep(1,N),Z)
  }
  else if(params$bias>1){
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
  V_offsets <- min( data$X[idx_not_missing] )
  data$X[idx_not_missing]<-data$X[idx_not_missing]-V_offsets+1
  # eventually, apply external transform
  #if(t%in%params){
    #work in logarithm space
    #data.X(:,r) = params.t_1{r}(data.X(:,r)); # work in logarithm space better
    #data.C(r) = params.ext_dataType{r};
    # ---To be completed ---
  #}
 
  func_bit<-rep(1,dim(data$X)[2])
  print(params)
  readline("Press return to continue")
  # call .Rcpp wrapper function
  setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/")
  library(RcppGSLExample)
 hidden<-IBPsampler(data$X,data$C,Z,unlist(params$bias),func_bit,unlist(params$s2u),unlist(params$s2B),unlist(params$alpha),unlist(params$Nsim),unlist(params$maxK), unlist(params$missing))
                    
  # From the posterior
 R<-rep(1,D)
  return(hidden)
}
 
  
  
  
  