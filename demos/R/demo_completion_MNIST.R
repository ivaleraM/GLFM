#demo_completion_MNIST.R
#'@description Demo for the MNIST dataset, data completion
#'@param alpha is the concentration parameter for the IBP
#'@param N is the number of datapoints
#'@param s2x is the noise variance
#'@param Z: NxK matrix of feature patterns
#'@param Niter: number of iterations for the Gibbs sampler
#'@param maxR: maximum number of categories across all dimensions
#'@param transf_dummie is a 0-1 variable that indicates if there are data transformations
#'@param th is the threshold to filter out the variables that are not significant
#'@return Xcompl is the completed matrix with maximum a posteriori estimates where the
#'missing values were

demo_data_exploration_MNIST<-function(){
rm(list=ls())
graphics.off()
#setwd("src/GLFMR")
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
require(ggplot2)
source("GLFM_infer.R")
source("GLFM_complete.R")
source("GLFM_computePDF.R")
source("GLFM_plotPatterns.R")
source("GLFM_computeMAP.R")
datos_mninst <- read.csv(file="../../datasets/mnist_train_small_100.csv", header=TRUE,stringsAsFactors = FALSE, sep=",")
n_labels<-datos_mninst[,1]
Xauxi<-as.matrix(datos_mninst[,2:785],ncol=99,nrow= 784, byrow=TRUE)
Xfull<-matrix((Xauxi),nrow=784,ncol=99)
Xfull<-t(Xfull)
Xfull<-Xfull+1
perc_missing <- 0.3
missing_val <- -100
N<-dim(Xfull)[1]
D<-dim(Xfull)[2]
Cfull<-rep('n',D)
rand_mat<-matrix(runif(N*D), ncol=D, nrow=N)
mask_missings <- (rand_mat < perc_missing)
Xmiss <- Xfull
Xmiss[mask_missings] <- missing_val
data_mnist<-list("X"=Xmiss,"C"=Cfull)
param_names<-c("missing","s2B","alpha","Niter","maxK","bias","transf_dummie")
bias <- 0 
s2B <- 0.5      
alpha <- 10   
Niter <- 1   
maxK <- 100  
missing <- -100
transf_dummie<-FALSE
params<-list(missing,s2B,alpha,Niter,maxK,bias,transf_dummie)
names(params)<-param_names
output2<-GLFM_complete(data_mnist,list(c(),params))
data_mnist$C<-output2$data$C
# Visualization of a random image
idxs_rndrow<-sample(1:dim(Xfull)[1],28,replace = TRUE, prob=rep(1/dim(Xfull)[1],dim(Xfull)[1]))
idxs_rndcol<-sample(1:dim(Xfull)[2],28,replace = TRUE, prob=rep(1/dim(Xfull)[2],dim(Xfull)[2]))
pixels<-Xfull[idxs_rndrow,idxs_rndcol]
image(pixels,col = grey(seq(0, 1, length = 256)))
# Example with missing 
pixels<-Xmiss[idxs_rndrow,idxs_rndcol]
dev.new()
image(pixels,col = grey(seq(0, 1, length = 256)))
# Example with the completed matrix
pixels<-output2$X_compl[idxs_rndrow,idxs_rndcol]
dev.new()
image(pixels,col = grey(seq(0, 1, length = 256)))
return(list("Xcompl"=output2$X_compl))
}
