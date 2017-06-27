#' @param theta: 1*R vector
#' @param Zp: 1*K vector
#' @param B: K*D matrix
pdf_o<-function(Zp,B,theta,s2Y){
  R<-length(theta)+1
  pdf<-rep(0,R)
  for(r in 1:R){
    if(r==1){
      a <- pnorm(theta[r], Zp%*%B, sqrt(s2Y) )
      b <- 0
    }
    else if(r==R){
      a <- 1
      b <-pnorm(theta[r-1], Zp%*%B, sqrt(s2Y) )
    }
    else{
      a <- pnorm(theta[r], Zp%*%B, sqrt(s2Y) )
      b <-pnorm(theta[r-1], Zp%*%B, sqrt(s2Y) )
    }
    pdf[r] = a-b
  }
}

