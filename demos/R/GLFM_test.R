#' @description Tiny example to check that the package with the Ccode is properly installed

rm(list=ls())
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("GLFM_complete.R")
source("GLFM_computePDF.R")
source("aux/init_default_params.R")
print(" Script to check call to GLFMR library")
D<-5
N <-3
X<-matrix(c(c(1.0, 1, -0.3, 1, 1),c(6.3, 2, 3.8, 23, 1),c(11, 3, 4.1, 4, 2)),nrow =D,ncol=N, byrow=TRUE)
C<-c('p','o','G','N','c')
data_test<-list("X"=t(X),"C"=C)
K <-2
Z<-matrix(c(c(1.0,0),c(1,1),c(1,1)),nrow=N,ncol=K)
print(list("The data matrix is:"=X,"The hidden matrix Z is:"=Z))
print("Checks if we can call C correctly:")
output<-GLFM_infer(data_test,list(Z,c()))
print(list("The inferred Z matrix is:"=output$hidden$Z,"the B matrix"=output$hidden$B))
print("Succesful")
