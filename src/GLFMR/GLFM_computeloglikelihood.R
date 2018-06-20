#'@description This function calculates the log-lik of
#' the held-out data in the GLFM. The function uses the trained
#' parameters,
#' to calculate the PDF at the points of missing values, using the true values.
#' Inputs:
#'@param data is a list with X and C
#'@param X: N*D data matrix X
#'@param C: 1xD string with data types, D = number of dimensions

GLFM_computeloglikelihood<-function(data,Zp,hidden,params,varargin){
  
  source("mapping_functions/f_o.R")
  source("mapping_functions/f_g.R")
  source("mapping_functions/f_c.R")
  source("mapping_functions/f_p.R")
  source("mapping_functions/f_o.R")
  source("mapping_functions/f_n.R")
  
  if (length(varargin) == 1){
    idxsD <- varargin[1]
  }
  else if(length(varargin) > 1){
    stop('Too many input arguments')
  }
  else{
    idxsD <- 1:length(hidden$B)
  }
  P <- dim(Zp)[1]
  K <-dim(hidden$B[[1]])[1]
  #readline("press return to continue")
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  D <- dim(data$X)[2]
  N <- dim(data$X)[1]
  lik<-matrix(0,N,D)
  for(d in 1:D){
    
  }
  
  
  
  
}