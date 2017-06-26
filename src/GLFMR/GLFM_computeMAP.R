#' Function to generate the MAP solution corresponding to patterns in Zp
#' Inputs:
#'@param C: 1xD string with data types, D = number of dimensions
#'@param Zp: PxK matrix of feature patterns for which it computes the MAP estimate
#'    (P is the number of feature patterns)
#' @param hidden: structure with latent variables learned by the model
#'@param B: latent feature list with D elements of size ( K * maxR)  where
#'@param D: number of dimensions
#'@param K: number of latent variables
#'@param maxR: maximum number of categories across all dimensions
#'@param  mu: 1*D shift parameter
#' @param w:  1*D scale parameter
#' @param  theta: D*maxR matrix of auxiliary vars (for ordinal variables)
  # ----------------(optional) ------------------
  #' @param     - idxsD: dimensions to infer
#' @return X_map: P*Di matrix with MAP estimate where Di = length(idxsD)
 
#It does not allow external transformations yet!

GLFM_computeMAP<-function(C,Zp,hidden,params,varargin){
  # We need to call the transformation functions! 
  source("f_o.R")
  source("f_g.R")
  source("f_c.R")
  source("f_p.R")
  source("f_o.R")
  source("f_n.R")
  if (length(varargin) == 1){
    idxsD <- varargin[1]
  }
  else if(length(varargin) > 1){
    stop('Too many input arguments')
  }
  else{
    idxsD <- 1:dim(hidden$B)[1]
  }
  P <- dim(Zp)[1]
  K <- dim(hidden$B)[2]
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  X_map<-matrix(0,nrow=P,ncol=length(idxsD))
  # For each dimension
  for(dd in 1:length(idxsD)){ # for each dimension
  d <- idxsD(dd)
  #('g','p','n','c','o')
  
  switch(C(d),'g'={X_map[,dd]<-f_g(Zp*unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])},
         'p'={X_map[,dd]<-f_p(Zp%*%unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])},
         'n'={X_map[,dd]<-f_n(Zp%*%unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])},
         'o'={X_map[,dd]<-f_o(Zp%*%unlist(hidden$B[d]),hidden$theta[d,1:(hidden$R[d]-1)])},
         'c'={X_map[,dd]<-f_c(Zp%*%unlist(hidden$B[d]))},
         stop('Unknown data type'))
  if(sum(is.nan(X_map[,dd])) > 0){
    warning('Some values are nan!') 
  }
    
  }
}


  
   