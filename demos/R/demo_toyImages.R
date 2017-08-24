# Demo toyImages
# image(m) es el analogo a imagesc en matlab

rm(list=ls())
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("generate_toy_images.R")
source("init_default_params.R")
# Generative model
N<-1000
s2x<-1

#Initialisation
data_gen<- generate_toy_images(N,s2x)
data<-list("X"=data_gen$X,"C"=data_gen$C)
#gT<-list(""=)
m0<-matrix(0,nrow=N,ncol=1)
Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.2,0.8)))
# define params
alpha <- 2   # Concentration parameter of the IBP
Niter <- 100  # Number of iterations for the Gibbs sampler
maxK <- 10
transf_dummie <-FALSE
params<-list("alpha"=alpha,"Niter"=Niter,"maxK"=maxK,"transf_dummie"=transf_dummie)
# Inference
output <- GLFM_infer(data, list(Z,params))
Kest <-dim(output$hidden$B[[1]])[1]
Zp <-diag(Kest)
X_map <- GLFM_computeMAP(data$C, Zp, output$hidden, output$params,c())
# Plot, ground truths
for(k in 1:dim(data_gen$gTB)[1]){
image(matrix(data_gen$gTB[k,],nrow=6,ncol=6))
  readline("Press return to continue")
}
# Plot inferred 
#X_map[which(X_map)<0]<-0
for(k in 1:Kest){
  image(matrix(X_map[k,],nrow=6,ncol=6))
  readline("Press return to continue")
}

