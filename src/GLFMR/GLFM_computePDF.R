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
  P <- dim(Zp)[1]
  K <-dim(hidden$B[[1]])[1]
  #readline("press return to continue")
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  if(data$C[d] == 'g' || data$C[d] == 'p'){
    if((params$numS %in% params)){
      numS<-params$numS
      xd <- seq(mm,MM,length.out=params$numS)
      print("xd is a linspace")
    }
    else{
      numS <-100
      xd <- seq(mm,MM,length.out=numS)
    }
  }
  else if(data$C[d] == 'n'){
    xd <-mm:MM
    numS <-length(xd)
  }
  else{
    xd <- unique(XXd[idxs_nonnans])
    print("xd is unique values")
    numS<-length(xd)
  #  params<-append("numS"=length(xd),params)
  }
  pdf_val <-matrix(0,P,numS)
  #print(list(numS,pdf_o(Zp[1,],hidden$B[[2]],hidden$theta[2,1:(2)],hidden$s2y[2])))
  # 2-8 equations in the paper for pdfs
  for(p in 1:P){
    switch(data$C[d],'g'={pdf_val[p,]<-pdf_g(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d],params)},
           'p'={pdf_val[p,]<-pdf_p(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d],params)},
           'n'={pdf_val[p,]<-pdf_n(xd,Zp[p,],hidden$B[[d]],hidden$mu[d],hidden$w[d],hidden$s2y[d],params)},
           'c'={pdf_val[p,]<-pdf_c(Zp[p,],hidden$B[[d]],hidden$s2y[d])},
           'o'={pdf_val[p,]<-pdf_o(Zp[p,],hidden$B[[d]],hidden$theta[d,1:(hidden$R[d]-1)],hidden$s2y[d])},
           stop('Unknown data type'))
  }
  if(sum(is.nan(pdf_val)) > 0){
    stop('Some values are nan!')
  }
  if("transf_dummie" %in% names(params)){
    if(params$transf_dummie){
      xd <- params$t_inv(xd)
      pdf_val <- pdf_val*abs(params$dt_1(xd))
    }
  }
  return(list("pdf"=pdf_val,"xd"=xd))
  }