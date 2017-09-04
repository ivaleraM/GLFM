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
#' @return X_map: P*D matrix with MAP estimate where Di = length(idxsD)
 
#It does not allow external transformations yet!
# It is not reading the dimension of B properly, 5 element list where 
# hidden$B is a list where each element is a vector of size 2
GLFM_computeMAP<-function(C,Zp,hidden,params,varargin){
  # We need to call the transformation functions! 
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
  #print(dim(Zp)[2])
  K <-dim(hidden$B[[1]])[1]
  #readline("press return to continue")
  #print(Zp)
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  X_map<-matrix(0,nrow=P,ncol=length(idxsD))
  
  for(jj in 1:length(idxsD)){# for each dimension
   dd<-idxsD[jj]
  switch(C[dd],'g'={X_map[,jj] <- f_g(Zp%*%hidden$B[[dd]],hidden$mu[dd],hidden$w[dd])},
         'p'={X_map[,jj] <- f_p(Zp%*%hidden$B[[dd]],hidden$mu[dd],hidden$w[dd])},
         'n'={X_map[,jj] <- f_n(Zp%*%hidden$B[[dd]],hidden$mu[dd],hidden$w[dd])},
         'o'={X_map[,jj] <- f_o(Zp%*%hidden$B[[dd]],hidden$theta[dd,1:(hidden$R[dd]-1)])},
         'c'={X_map[,jj] <- f_c(Zp%*%hidden$B[[dd]])},
         stop('Unknown data type'))
    if(sum(is.nan(X_map[,jj])) > 0){
    warning('Some values are nan!') 
    }
  }
   if("transf_dummie" %in% names(params)){
    if(params$transf_dummie){
      if(is.list(params$t_1)==FALSE){
      X_map[,params$idx_transform] <-params$t_inv( X_map[,params$idx_transform])
      }else{
      for(ell in 1:length(params$t_inv)){
        X_map[,params$idx_transform[[ell]]] <-params$t_inv[[ell]]( X_map[,params$idx_transform[[ell]]])
        }
    }
    }
   }
  return(X_map)
  }



  
   