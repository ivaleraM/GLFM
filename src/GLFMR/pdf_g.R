pdf_g<-function(x,Zp,B,mu,w,s2Y,params){
  
 df_1<-function(x){
   w*(x-mu)
 }
 pdf <-pnorm(df_1(x),Zp * B,sqrt(s2Y + params$s2u))%*%w
}



