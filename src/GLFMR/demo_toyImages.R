# Demo toyImages
rm(list=ls())
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("generate_toy_images.R")

# Generative model
N<-1000
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
hidden <- GLFM_infer(data, list(Z,params))
Kest<-dim(hidden$B)[2]
Z_p <-diag(Kest)
X_map <- GLFM_computeMAP(data$C, Zp, hidden, params)
# Plot