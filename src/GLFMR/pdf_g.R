pdf_g<-function(xx,Zp,B,mu,w,s2Y,params){
source("df_1.R")
 pdf <-dnorm(df_1(xx,w,mu),Zp%*%B,sqrt(s2Y + params$s2u))
 #print(list(df_1(xx,w,mu),Zp%*%B,sqrt(s2Y + params$s2u),w))
  pdf<-pdf*w
  
 #pdf<-dnorm(df_1(xx,w,mu))
}



