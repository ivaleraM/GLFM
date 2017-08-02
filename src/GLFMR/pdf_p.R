#' 
pdf_p<-function(x,Zp, B, mu, w, s2Y){
source("f_p_1.R") 
  pdf_out <- 1/(2*pi*sqrt(s2Y + params$s2u))*exp( -1/(2*(s2Y + params$s2u))*(f_p_1(x, mu, w) - Zp %*%B)*as.vector(f_p_1(x, mu, w) - Zp %*%B))
  pdf_out<-pdf_out*as.vector(abs(df_p_1(x, mu, w)) )
 # print(list(pdf_out,abs(df_p_1(x, mu, w))))
  return(pdf_out)
}

