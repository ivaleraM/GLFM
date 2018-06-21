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
  Z_p <-hidden$Z
  # Deals with missing values
  X_aux<-data$X
  aa<-max(X_aux)
  X_aux[idx_missing] <- aa+1
  V_offset<-colMins(X_aux)
  V_offset_mat<-matrix(V_offset,nrow=N,ncol=D,byrow=TRUE)
  X_aux<-data$X-V_offset_mat+1
  idx_catord<-which(data$C=='c' | data$C=='o')
  if(length(idx_catord)>0){
    data$X[,idx_catord] <-X_aux[,idx_catord]
    bu<-apply(X_aux[,idx_catord,drop=FALSE], 2, function(x)length(unique(x)))
    idx_dat<-which(colMaxs(X_aux[,idx_catord,drop=FALSE])!=bu)
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
  
# if there is an external transformation change type of dimension d by external data type
    if( "transf_dummie" %in% names(params)){
      if(params2$transf_dummie){
        if(is.list(params2$t_1)==FALSE){
          data$X[,params2$idx_transform]<-params2$t_1(data$X[,params2$idx_transform])
          data$C[params2$idx_transform] <-params2$ext_datatype
        }else{
          for(ell in 1:length(params2$t_1)){
            data$X[,params2$idx_transform[[ell]]]<-params2$t_1[[ell]](data$X[,params2$idx_transform[[ell]]])
            data$C[params2$idx_transform[[ell]]] <-params2$ext_datatype[[ell]]
          }
        }
      }
    }
    
      #Find coordinates of missing values (NaN's are considered as missing)
      idxs_nonnans<-which(!is.nan(X_true), arr.in=TRUE)
      
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