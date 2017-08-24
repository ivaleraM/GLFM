#' @param theta: 1*R vector
#' @param Zp: 1*K vector
#' @param B: K*D matrix
pdf_o<-function(Zp,B,theta,s2Y){
  R<-length(theta)+1
  pdf_out<-rep(0,R)
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
    if(length(a)>1 ||length(b)>1){
    stop("a,b values should be scalar!")
    }
    pdf_out[r] = a-b
  }
  return(pdf_out)
}

