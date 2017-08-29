#' @description Function to complete a matrix that has missing values
#' @param data a list with the N*D matrix X and C 1*D vector
#' @param hidden feature assignment N*K matrix Z
#' @param params a list with the simulation parameters and hyperparameters
#' @return Xcompl the completed data matrix 
GLFM_complete<-function(data,varargin){
  source("aux/init_default_params.R")
  X_compl<-data$X
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
  if(sum(is.nan(data$X))==0 && (sum(data$X == params$missing) == 0 || is.na(sum(data$X == params$missing)))){
    print('The input matrix X has no missing values to complete.')
    Xcompl <- c()
  }
  else{
    print("press enter to continue")
    output <- GLFM_infer(data, list(Z,params))
    # NaN's are considered as missing
    data$X[which(is.nan(data$X))]=params$missing
    idxs_missing<-which(data$X == params$missing)
    print(idxs_missing)
    readline("press return to continue")
    # For each missing value, compute a MAP estimate
    for(j in 1:length(idxs_missing)){
    X_compl[idxs_missing[j]]<-GLFM_computeMAP(data$C, Z[idxs_missing[j]], output$hidden, output$params,idxs_missing[j])
    }
  }
  return(list("X_compl"=X_compl,"hidden"=hidden))
        
}



  
  
  
  
