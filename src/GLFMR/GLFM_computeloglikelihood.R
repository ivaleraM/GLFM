#'@description This function calculates the log-lik of
#' the held-out data in the GLFM. The function uses the trained
#' parameters,
#' to calculate the PDF at the points of missing values, using the true values.
#' Inputs:
#'@param data is a list with X and C
#'@param X: N*D data matrix X
#'@param C: 1xD string with data types, D = number of dimensions

GLFM_computeloglikelihood<-function(data,hidden,params,varargin){
  
  source("mapping_functions/f_o.R")
  source("mapping_functions/f_g.R")
  source("mapping_functions/f_c.R")
  source("mapping_functions/f_p.R")
  source("mapping_functions/f_o.R")
  source("mapping_functions/f_n.R")
  
  # Not sure if I need extra stuff in varargin yet
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
  # Copy of the data
  X_true <- data$X
  X_transformed <- data$X
  Z_p <-hidden$Z 
  for(d in 1:D){
    # if there is an external transformation change type of dimension d by external data type
    #
    #Find coordinates of missing values (NaN's are considered as missing)
    idxs_nonnans<-which(!is.nan(data$X), arr.in=TRUE)
    if(data$C[d] == 'g' || data$C[d] == 'p'){
      if((params$numS %in% params)){
        numS<-params$numS
        xd <- seq(mm,MM,length.out=params$numS)
        print(list("xd is a linspace, positive or real data"))
      }
      else{
        numS <-100
        xd <- seq(mm,MM,length.out=numS)
      }
    }
    else if(data$C[d] == 'n'){
      xd <-mm:MM
      numS <-length(xd)
      print("xd is a grid, count data")
    }
    else{
      #***difficult bit
      xd <- unique(X[idxs_nans[1,d])
      numS<-length(xd)
      print(list("xd are unique values, categorical or ordinal data"))
    }
  }
  # Gives the number of non-missing entries:
 rowsnum<-dim(idxs_nonnans)[1]
 lik<-rep(0,rowsnum)
  for(ell in rowsnum){
    n_idx = idxs_nonnans[ell,][1]
    d_idx = idxs_nonnans[ell,][2]
    xd = X_true[n_idx,d_idx]
    switch(data$C[d],'g'={lik[ell]<-pdf_g(xd,Zp[n_idx,],hidden$B[[d_idx]],hidden$mu[d_idx],hidden$w[d_idx],hidden$s2y[d_idx],params)},
           'p'={lik[ell]<-pdf_p(xd,Zp[n_idx,],hidden$B[[d_idx]],hidden$mu[d_idx],hidden$w[d],hidden$s2y[d_idx])},
           'n'={lik[ell]<-pdf_n(xd,Zp[n_idx,],hidden$B[[d_idx]],hidden$mu[d_idx],hidden$w[d],hidden$s2y[d_idx],params)},
           'c'={lik[ell]<-pdf_c(Zp[n_idx,],hidden$B[[d_idx]],hidden$s2y[d_idx])},
           'o'={lik[ell]<-pdf_o(Zp[n_idx,],hidden$B[[d_idx]],hidden$theta[d,1:(hidden$R[d_idx]-1)],hidden$s2y[d_idx])},
           stop('Unknown data type'))
    if(sum(is.nan(lik[ell])) > 0){
      print(data$C[d])
      stop('Some values are nan!')
    }
  }
 
  
  
  
}