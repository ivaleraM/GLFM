pdf_g<-function(xx,Zp,B,mu,w,s2Y,params){
source("mapping_functions/df_1.R")
 pdf <-dnorm(df_1(xx,w,mu),Zp%*%B,sqrt(s2Y + params$s2u))
 print(dnorm(df_1(xx,w,mu),Zp%*%B,sqrt(s2Y + params$s2u)))
  pdf<-pdf*w
  
 #pdf<-dnorm(df_1(xx,w,mu))
}



