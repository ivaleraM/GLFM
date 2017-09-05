# Another tiny example

rm(list=ls())
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("GLFM_complete.R")
source("GLFM_computePDF.R")
source("aux/init_default_params.R")
 X <- matrix(rnorm(10,0,1),nrow=10,ncol=5)
 C <- c('g','g','g','g','g')
 data<-list("X"=X,"C"=C)
 output<-GLFM_infer(data,c())
 Kest<-dim(output$hidden$B[[1]])[1]
 Zp <-diag(Kest)
 X_map <- GLFM_computeMAP(data$C, output$hidden$Z, output$hidden, output$params,c())
 D<-dim(X)[2]
 for(dd in 1:D){
  pdf_val<-GLFM_computePDF(data,Zp,output$hidden,output$params,dd) 
  print(pdf_val)
  readline("Press return to continue")
 }
 X[1,5]<--1
 X[10,5]<--1
 data<-list("X"=X,"C"=C)
output2<-GLFM_complete(data,list(c(),c()))
 