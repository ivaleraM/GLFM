#' Inputs:
#'@param data is a list with X and C
#'@param X: N*D data matrix X
#'@param C: 1xD string with data types, D = number of dimensions
#'@param Zp: P*K matrix of patterns (P is the number of patterns)
#'@param hidden: list with the following latent variables learned by the model (B,Z):
#'@param B: latent feature list with D matrices of size  K * maxR  where
#'@param D: number of dimensions
#'@param K: number of latent variables
#'@param maxR: maximum number of categories across all dimensions
#'@param Z: PxK matrix of feature patterns
#'@param params: list with parameters (mu,w,theta)
#'@param  mu: 1*D shift parameter
#'@param w:  1*D scale parameter
#'@param  theta: D*maxR matrix of auxiliary vars (for ordinal variables)
#'Outputs:
#'@param d: dimension for which the pdf is computed
#'@return A list with two vectors (xd,pdf_val) where
#'@return xd 1*NumS, NumS points are to be considered
#'@return pdf P*numS matrix where P is the number of patterns 

GLFM_computePDF<-function(data,Zp,hidden,params,d){
  source("pdf_functions/pdf_g.R")
  source("pdf_p.R")
  source("pdf_functions/pdf_n.R")
  source("pdf_functions/pdf_c.R")
  source("pdf_functions/pdf_o.R")
  source("df_p_1.R")
  XXd<-data$X[,d]
  idxs_nans <- which(is.nan(XXd))
  if(length(idxs_nans) > 0){
    XXd[idxs_nans] = params$missing
  }
  idxs_nonnans<-setdiff(1:(length(XXd)),idxs_nans)
  mm <- min(XXd[idxs_nonnans])
  MM <- max(XXd[idxs_nonnans]) 
  # External transformation case
  if("transf_dummie" %in% names(params) ){
    if(is.list(params$t_1)==FALSE){
    if(params$transf_dummie && d %in% params$idx_transform){
      mm <- params$t_1(mm)
      MM <- params$t_1(MM)
      }
    }else{
      for(ell in 1:length(params$t_1)){
      if(d %in% params$idx_transform[[ell]]){
      mm <- params$t_1[[ell]](mm)
      MM <- params$t_1[[ell]](MM)
          }
        }
      }
    }
  P <- dim(Zp)[1]
  K <-dim(hidden$B[[1]])[1]
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
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
    xd <- unique(XXd[idxs_nonnans])
    print(list("xd are unique values, categorical or ordinal data"))
    
    numS<-length(xd)
  }
  pdf_val <-matrix(0,P,numS)
  
  for(p in 1:P){
    switch(data$C[d],'g'={pdf_val[p,]<-pdf_g(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d],params)},
           'p'={pdf_val[p,]<-pdf_p(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d])},
           'n'={pdf_val[p,]<-pdf_n(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d],params)},
           'c'={pdf_val[p,]<-pdf_c(Zp[p,],hidden$B[[d]],hidden$s2y[d])},
           'o'={pdf_val[p,]<-pdf_o(Zp[p,],hidden$B[[d]],hidden$theta[d,1:(hidden$R[d]-1)],hidden$s2y[d])},
           stop('Unknown data type'))
  }
  if(sum(is.nan(pdf_val)) > 0){
    print(data$C[d])
    stop('Some values are nan!')
  }
  if("transf_dummie" %in% names(params)){
    if(is.list(params$t_1)==FALSE){
   if(params$transf_dummie && d %in% params$idx_transform){
      xd <- params$t_inv(xd)
      pdf_val<-pdf_val%*%diag(abs(params$dt_1(xd)))
    }
    }else{
      for(ell in 1:length(params$t_1)){
        if(d %in% params$idx_transform[[ell]]){
          xd <- params$t_inv[[ell]](xd)
          pdf_val<-pdf_val%*%diag(abs(params$dt_1[[ell]](xd)))
        }
      }
  }
  }
  return(list("pdf"=pdf_val,"xd"=xd))
  }