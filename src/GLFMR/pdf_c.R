#' @param B: K*R matrix
#' @param Zp: 1*K, where K: number of latent features
#' @return the pdf for categorical data
 pdf_c<-function(Zp, B, s2Y){
  numMC <- 1000 #% number of Monte Carlo samples to approximate the Expectation
  RR <- dim(B)[2]
  pdf_out <-rep(0,RR)
  u <- sqrt(s2Y) * rnorm(numMC)
  for(rr in 1:RR){
    aux<-rep(1,numMC)
    for(ell in 1:RR){
      if(ell==rr){
        next()}
        #aux <- aux%*%pnorm( u + Zp*(B[,rr] - B[,ell] ))
        aux <- aux *as.vector(pnorm( u + Zp%*%(B[,rr] - B[,ell] )))
      # print(Zp%*%(B[,rr] - B[,ell] ))
      }
    pdf_out[rr]<-mean(aux)
    #print(mean(aux))
  }
  pdf_out<-pdf_out/sum(pdf_out)
  return(pdf_out)
 }
