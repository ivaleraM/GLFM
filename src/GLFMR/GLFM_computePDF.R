#' @param data: is a list with N*D matrix X and 1*D vector C
#' @param Zp: P*K matrix of patterns (P is the number of patterns)
#' @param hidden$B: is a list of D elements: each element is a K*maxR vector
#' @param hidden$mu: 1*D shift parameter
#' @param hidden$w: 1*D scale parameter 
#' @param hidden$theta: D*maxR matrix of auxiliary variables for ordinal data
#' where maxR is the maximum number of categories
#' @return xd 1*NumS, NumS points are to be considered
#' @return pdf P*numS matrix where P is the number of patterns 

GLFM_computePDF<-function(data,Zp,hidden,params,d){
  idxs_nans <- which(is.nan(data$X[,d]))
  if(length(idxs_nans) > 0){
    data$X[idxs_nans,d] = params$missing
  }
  idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
  mm <- min(data$X[idxs_nonnans,d])
  M <- max(data$X[idxs_nonnans,d])
  # Add external transformation case
  B_aux<-matrix(unlist(hidden$B),nrow=dim(hidden$B)[1],ncol=dim(Zp)[2],byrow=TRUE)
  P <- dim(Zp)[1]
  K <-dim(B_aux)[2]
  #readline("press return to continue")
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  if(data$C[d] == 'g' || data$C[d] == 'p'){
    if((param$numS %in% params)==FALSE){
      params<-append("numS"=100,params)
      xd <- seq(mm,MM,length=numS)
    }
  }
  else if(data$C[d] == 'n'){
    xd <-mm:MM
    params<-append("numS"=length(xd),params)
  }
  else{
    xd <- unique(data$X[idxs_nonnans,d])
    params<-append("numS"=length(xd),params)
  }
  pdf <-matrix(0,P,params$numS)
  for(p in 1:P){
    switch(data$C[d],'g'={pdf[p,]<-pdf_g(Zp[p,]%*%B_aux[d,],hidden$mu[d],hidden$w[d],params$s2y,params)},
           'p'={pdf[p,]<-pdf_p(Zp[p,]%*%B_aux[d,],hidden$mu[d],hidden$w[d],params$s2y,params)},
           'n'={pdf[p,]<-pdf_n(Zp[p,]%*%B_aux[d,],hidden$mu[d],hidden$w[d],params$s2y,params)},
           'o'={pdf[p,]<-pdf_c(Zp[p,]%*%B_aux[d,1:(hidden$R[d]-1)],params$s2y)},
           'c'={pdf[p,]<-pdf_o(Zp[p,]%*%B_aux[d,],hidden$theta[d,1:(hidden$R[d]-1)],params$s2y,params)},
           stop('Unknown data type'))
  }
  if(sum(is.nan(pdf)) > 0){
    stop('Some values are nan!')
  }
  # External transformation case is missing
  }