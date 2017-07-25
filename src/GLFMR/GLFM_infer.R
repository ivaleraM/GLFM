#' @param data is a list with X and C
#' @param varargin is a list of lists containing hidden variables and parameters 
#' @return A list with posterior draws

GLFM_infer<-function(data,varargin){
  require(matrixStats)
 
  #print(length(varargin[[1]]))
  varargin_size<-length(varargin)
  if(varargin_size==0){
    Z<-c()
    params2 <- init_default_params(data, c())
  }else if(varargin_size==1){
      Z<-unlist(varargin[[1]])
      params2 <- init_default_params(data, c())
  }else if(varargin_size==2){
    Z<-unlist(varargin[[1]])
    params2 <- init_default_params(data, varargin[[2]])
  }else{
    stop("Incorrect number of input parameters: should be 0, 1 or 2")
  }
  D <- dim(data$X)[2]
  N <- dim(data$X)[1]
 
  if(length(Z)==0){
    m0<-matrix(0,N,2)
    Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
  }
  if(params2$bias == 1 && length(params2$bias)>0){
    Z <-cbind(rep(1,N),Z)
  }
  else if(params2$bias>1  && length(params2$bias)>0){
    stop("There is more than 1 bias specified, but the structure Z has not been initialized") 
  }
  # replace missing values
  idx_missing<-which(is.nan(data$X))
  data$X[idx_missing] <- params2$missing
  # Change labels of categorical or ordinal data such that the minimum value is 1
  idx_missing<-which((data$X)==params2$missing)
  X_aux<-data$X
  aa<-max(X_aux)
  X_aux[idx_missing] <- aa+1
  V_offset<-colMins(X_aux)
  V_offset_mat<-matrix(V_offset,nrow=N,ncol=D,byrow=TRUE)
  X_aux<-data$X-V_offset_mat+1
  idx_catord<-which(data$C=='c' | data$C=='o')
  if(length(idx_catord)>0){
    data$X[,idx_catord] <-X_aux[,idx_catord]
    bu<-apply(X_aux[,idx_catord], 2, function(x)length(unique(x)))
    #while(any(colMaxs(X_aux[,idx_catord])!=bu)){
      idx_dat<-which(colMaxs(X_aux[,idx_catord])!=bu)
      if(length(idx_dat)>0){
      for(ii in 1:length(idx_dat)){
      idxs_bad<-which(X_aux[,idx_dat[ii]]>bu[idx_dat[ii]])
      while(length(idxs_bad)>0){
        X_aux[idxs_bad,idx_dat[ii]]<-X_aux[idxs_bad,idx_dat[ii]]-1
        idxs_bad<-which(X_aux[,idx_dat[ii]]>bu[idx_dat[ii]])
          }
        }
      }
    data$X[,idx_catord]<-X_aux[,idx_catord]
    data$X[idx_missing]<-params2$missing 
      }
  
  
  # eventually, apply external transform
  if( "transf_dummie" %in% names(params2)){
    if(params2$transf_dummie){
     data$X[,params2$idx_transform]<-params2$t_1(data$X[,params2$idx_transform])
  #  data$C[params2$idx_transform] <-params2$ext_datatype
    }
  }
  
  func_bit<-rep(1,dim(data$X)[2])
  setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/")
  library(RcppGSLExample)
  #print(params2)
  print(list(length(data$C),dim(data$X)))
  readline("Press return to continue")
  # call .Rcpp wrapper function
  
  hidden<-IBPsampler(t(data$X),(data$C),t(Z),params2$bias,func_bit,params2$s2u,params2$s2B,
                    params2$alpha,params2$Niter,params2$maxK,params2$missing)
 #print(hidden)
 #readline("press return to continue")
  R<-rep(1,D)
  if(length(idx_catord)>0){
    X_aux<-data$X
    aa<-min(X_aux)
    X_aux[idx_missing] <- aa-1
    V_offset<--colMins(-X_aux)
    R[idx_catord]<-V_offset[idx_catord]
  }
  hidden$Z<-t(hidden$Z)
 hidden<-append(hidden, list("R"=R))
 setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
  return(list("hidden"=hidden,"params"=params2))
}
 
  
  
  
  