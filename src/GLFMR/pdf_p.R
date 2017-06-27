#' 
pdf_p<-function(x,Zp, B, mu, w, s2Y, params){
source("f_p_1.R") 
  pdf <- 1/(2*pi*sqrt(s2Y + params$s2u))*exp( -1/(2*(s2Y + params$s2u))*(f_p_1(x, mu, w) - Zp %*%B)^2 ) 
  pdf<-pdf* abs(f_p_1(x, mu, w)) 
}

