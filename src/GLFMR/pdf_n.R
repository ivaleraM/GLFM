#' @param B: K*R matrix
#' @param Zp: 1*K, where K: number of latent features
#' @return the pdf for nominal data
pdf_n<-function(x,Zp, B, mu, w, s2Y, params){
  source("f_p_1.R")
 pdf_out<-pnorm(f_p_1(x+1,mu,w),Zp%*%B,sqrt(s2Y))- pnorm(f_p_1(x,mu,w),Zp%*%B,sqrt(s2Y))
 return(pdf_out)
}

