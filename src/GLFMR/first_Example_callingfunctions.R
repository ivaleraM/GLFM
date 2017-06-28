# Tiny example
 X <- matrix(rnorm(10,0,1),nrow=10,ncol=5)
 C <- c('g','g','g','g','g')
 m0 <- matrix(0,10,2)
 Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
 data<-list("X"=X,"C"=C)
 # GLFM_infer, c() es para un vector vacio:
 output<-GLFM_infer(data,c())
 # GLFM_computeMAP
 Kest<-dim(output$hidden$B)[1]
 Zp <-diag(Kest)
 X_map <- GLFM_computeMAP(data$C, Zp, output$hidden, output$params,c())
 #GLFM_complete
 X[1,5]<-NaN
 X[10,5]<-NaN
 data<-list("X"=X,"C"=C)
output2<-GLFM_complete(data,c())
 # Aqui tengo duda sobre como llamar a GLFM_computeMAP