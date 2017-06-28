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
  source("pdf_g.R")
  source("pdf_p.R")
  source("pdf_n.R")
  source("pdf_c.R")
  source("pdf_o.R")
  XXd<-data$X[,d]
  #print(XXd)
  idxs_nans <- which(is.nan(XXd))
  if(length(idxs_nans) > 0){
    XXd[idxs_nans] = params$missing
  }
  idxs_nonnans<-setdiff(1:(length(XXd)),idxs_nans)
  mm <- min(XXd[idxs_nonnans])
  MM <- max(XXd[idxs_nonnans])
  # Add external transformation case
  B_aux<-matrix(unlist(hidden$B),nrow=dim(hidden$B)[1],ncol=dim(Zp)[2],byrow=TRUE)
  P <- dim(Zp)[1]
  K <-dim(B_aux)[2]
  #readline("press return to continue")
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  if(data$C[d] == 'g' || data$C[d] == 'p'){
    if((params$numS %in% params)){
      xd <- seq(mm,MM,length.out=params$numS)
      print("xd is a linspace")
    }
    else{
      params<-append("numS"=100,params)
      xd <- seq(mm,MM,length.out=params$numS)
    }
  }
  else if(data$C[d] == 'n'){
    xd <-mm:MM
    params<-append("numS"=length(xd),params)
  }
  else{
    xd <- unique(XXd[idxs_nonnans])
    print("xd is unique values")
    params<-append("numS"=length(xd),params)
  }
  pdf_val <-matrix(0,P,params$numS)
  #print(xd)
  for(p in 1:P){
    switch(data$C[d],'g'={pdf_val[p,]<-pdf_g(xd,Zp[p,],B_aux[d,],hidden$mu[d],hidden$w[d],hidden$s2y,params)},
           'p'={pdf_val[p,]<-pdf_p(xd,Zp[p,],B_aux[d,],hidden$mu[d],hidden$w[d],hidden$s2y,params)},
           'n'={pdf_val[p,]<-pdf_n(xd,Zp[p,],B_aux[d,],hidden$mu[d],hidden$w[d],hidden$s2y,params)},
           'o'={pdf_val[p,]<-pdf_c(xd,Zp[p,],B_aux[d,1:(hidden$R[d]-1)],hidden$s2y)},
           'c'={pdf_val[p,]<-pdf_o(xd,Zp[p,],B_aux[d,],hidden$theta[d,1:(hidden$R[d]-1)],hidden$s2y,params)},
           stop('Unknown data type'))
  }
  if(sum(is.nan(pdf_val)) > 0){
    stop('Some values are nan!')
  }
  return(pdf_val)
  # External transformation case is missing
  }