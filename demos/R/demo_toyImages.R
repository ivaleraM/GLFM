#'@description Demo for the toyImages dataset
#'@param alpha is the concentration parameter for the IBP
#'@param N is the number of datapoints
#'@param s2x is the noise variance
#'@param Z: NxK matrix of feature patterns
#'@param Niter: number of iterations for the Gibbs sampler
#'@param maxR: maximum number of categories across all dimensions
#'@param transf_dummie is a 0-1 variable that indicates if there are data transformations
#'@return A list of two elements, output is a list with the output$hidden and output$params
#' lists and Xmap is the maximum a posteriori estimate

demo_toyImages<-function(){
rm(list=ls())
graphics.off()
#Duda de como cambiar a directorio relativo
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("generate_toy_images.R")
source("init_default_params.R")
N<-1000
s2x<-1
#Initialisation
data_gen<- generate_toy_images(N,s2x)
data<-list("X"=data_gen$X,"C"=data_gen$C)
m0<-matrix(0,nrow=N,ncol=1)
Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.2,0.8)))
alpha <- 2   
Niter <- 100  
maxK <- 10
transf_dummie <-FALSE
params<-list("alpha"=alpha,"Niter"=Niter,"maxK"=maxK,"transf_dummie"=transf_dummie)
# Inference
output <- GLFM_infer(data, list(Z,params))
Kest <-dim(output$hidden$B[[1]])[1]
Zp <-diag(Kest)
X_map <- GLFM_computeMAP(data$C, Zp, output$hidden, output$params,c())
for(k in 1:dim(data_gen$gTB)[1]){
image(matrix(data_gen$gTB[k,],nrow=6,ncol=6))
  readline("Press return to continue")
}
for(k in 1:Kest){
  image(matrix(X_map[k,],nrow=6,ncol=6))
  readline("Press return to continue")
}
return(list("output"=output,"Xmap"=Xmap))
}

