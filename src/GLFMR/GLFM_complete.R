#' @description Function to complete a matrix that has missing values
#' @param data a list with the N*D matrix X and C 1*D vector
#' @param hidden feature assignment N*K matrix Z
#' @param params a list with the simulation parameters and hyperparameters
#' @return Xcompl the completed data matrix 
GLFM_complete<-function(data,varargin){
  varargin_size<-length(varargin)
  if(varargin_size==0){
    Z<-c()
    params <- init_default_params(data, c())
    print("Case 0")
  }
  else if(varargin_size<3)
  {
    Z<-unlist(varargin[1])
    
    switch(varargin_size,params <- init_default_params(data, c()) ,params <- init_default_params(data, varargin[2]))
  }
  else
  {
    stop("Incorrect number of input parameters: should be 0, 1 or 2")
  }
  if(sum(is.nan(data$X)==0 || sum(data$X == params$missing) == 0)){
    prin('The input matrix X has no missing values to complete.')
    Xcompl = c()
  }
  else{
    output <- GLFM_infer(data, list(Z,params))
    # NaN's are considered as missing
    data$X[which(is.nan(data$X))]=params$missing
    idxs_missing<-which(data$X == params$missing)
    # For each missing value, compute a MAP estimate
    for(j in 1:length(idxs_missing)){
    X_aux<-GLFM_computeMAP(data$C, Z[idxs_missing[j]], output$hidden, output$params,idxs_missing[j])
    }
  }
        
}



  
  
  
  
