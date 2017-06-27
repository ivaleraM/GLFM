#' @param B: K*R matrix
#' @param Zp: 1*K, where K: number of latent features
#' @return the pdf for categorical data
 pdf_c<-function(Zp, B, s2Y){
  numMC <- 1000 #% number of Monte Carlo samples to approximate the Expectation
  R <- dim(B)[2]
  pdf <-rep(0,R)
  u <- sqrt(s2Y) * rnorm(numMC)
  for(rr in 1:R){
    aux<-rep(1,numMC)
    for(jj in 1:R){
      if(jj!=rr){
        aux <- aux%*%pnorm( u + Zp*(B[,rr] - B[,jj] ))
      }
    }
    pdf[rr]<-mean(aux)
  }
  pdf<-pdf/sum(pdf)
  return(pdf)
 }
