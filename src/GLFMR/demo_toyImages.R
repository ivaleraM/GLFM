# Demo toyImages
# image(m) es el analogo a imagesc en matlab

rm(list=ls())
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("generate_toy_images.R")

# Generative model
N<-200
s2x<-0.5

#Initialisation
data_gen<- generate_toy_images(N,s2x)
data<-list("X"=data_gen$X,"C"=data_gen$C)
#gT<-list(""=)
m0<-matrix(0,nrow=N,ncol=1)
Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
# define params
alpha <- 2   # Concentration parameter of the IBP
Niter <- 100  # Number of iterations for the Gibbs sampler
maxK <- 10
params<-list("alpha"=alpha,"Niter"=Niter,"maxK"=maxK)
# Inference
output <- GLFM_infer(data, list(Z,params))
Kest<-dim(output$hidden$B)[1]
Zp <-diag(Kest)
X_map <- GLFM_computeMAP(data$C, Zp, output$hidden, output$params,c())
# Plot, ground truths
for(k in 1:dim(data_gen$gTB)[1]){
image(matrix(data_gen$gTB[1,],nrow=6,ncol=6))
}
# Plot inferred 
for(k in 1:Kest){
  image(matrix(X_map[k,],nrow=6,ncol=6))
  readline("Press return to continue")
}

