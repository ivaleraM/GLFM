pdf_g<-function(x,Zp,B,mu,w,s2Y,params){
source("df_1.R")
 pdf <-pnorm(df_1(x,w,mu),Zp%*%B,sqrt(s2Y + params$s2u))%*%w
}



