#' @description Function to complete a data matrix X that has missing values
#' Inputs:
#'@param data is a list with the data matrix X and a string of data types C
#'@param X: N*D data matrix X
#'@param C: 1xD string with data types, D = number of dimensions
#' # ----------------(optional) ------------------
#' @param varargin is a list of lists containing 2 lists hidden and params 
#'@param hidden: list with the following latent variables (Z,B):
#'@param Z: NxK matrix of feature patterns
#'@param B: latent feature list with D matrices of size  K * maxR  where
#'@param D: number of dimensions
#'@param K: number of latent variables
#'@param maxR: maximum number of categories across all dimensions
#'@param params: list with parameters (mu,w,theta)
#'@param  mu: 1*D shift parameter
#'@param w: 1*D scale parameter
#'@param theta: D*maxR matrix of auxiliary vars (for ordinal variables)
#'@param s2y: 1*D per dim inferred noise variance for pseudo observations
#' @param hidden feature assignment N*K matrix Z
#' @param params a list with the simulation parameters and hyperparameters
#' Outputs:
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
  if(sum(is.nan(data$X))==0 && (sum(data$X == params2$missing) == 0 || is.na(sum(data$X == params2$missing)))){
    print('The input matrix X has no missing values to complete.')
    Xcompl <- c()
  }
  else{
    output <- GLFM_infer(data, list(Z,params2))
    # NaN's are considered as missing
    data$X[which(is.nan(data$X))]=params2$missing
    idxs_missing<-which(data$X == params2$missing, arr.in=TRUE)
    X_compl<-data$X
    for(j in 1:dim(idxs_missing)[1]){
       X_compl[idxs_missing[j,1],idxs_missing[j,2]]<-GLFM_computeMAP(data$C,output$hidden$Z[idxs_missing[j,1],,drop=FALSE], output$hidden, output$params,idxs_missing[j,2])
    }
  }
  return(list("X_compl"=X_compl,"hidden"=output$hidden))
        
}



  
  
  
  
