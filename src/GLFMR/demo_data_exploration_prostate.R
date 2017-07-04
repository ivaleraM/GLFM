# demo_data_exploration_prostate
rm(list=ls())
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/src/GLFMR")
require(R.matlab)
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
datos_prostate<-readMat('prostate_v3.mat')
source("init_default_params.R")
Xauxi <- as.matrix(unlist(datos_prostate$data[2,1,1]),ncol=16,nrow= 502, byrow=TRUE)
Xfull<-aux<-matrix(Xauxi,nrow=502,ncol=16)
C<-unlist(datos_prostate$data[3,1,1],use.names = FALSE)
Cfull<-strsplit(as.character(C), "")
cat_labels_full <-unlist(datos_prostate$data[4,1,1],use.names = FALSE)
y_labels_full<-unlist(datos_prostate$data[5,1,1],use.names = FALSE)
y_labels_long_full<-unlist(datos_prostate$data[6,1,1],use.names = FALSE)
idx_toKeep <- c(1, 2, 4,13, 15)
X<-Xfull[,idx_toKeep]
C<-Cfull[[1]][idx_toKeep]
cat_labels<-cat_labels_full[idx_toKeep]
y_labels<-y_labels_full[idx_toKeep]
y_labels_long<-y_labels_long_full[idx_toKeep]
#params
param_names<-c("missing","s2u","s2B","alpha","Niter","maxK","bias")
missing<--1
s2u<-0.005
s2B<-1
alpha<-1
Niter<-100
maxK<-10
bias<-1
params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias)
names(params)<-param_names
N<-dim(X)[1]
#m0<-matrix(0,N,2)
#Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
#if(params$bias == 1 && length(params$bias)>0){
#  Z <-cbind(rep(1,N),Z)
#}
# Inference
Z<-c()
data<-list("X"=X,"C"=C)
output <- GLFM_infer(data, list(Z,params))
#Predict MAP estimate for the whole matrix X
X_map <- GLFM_computeMAP(data$C, output$hidden$Z, output$hidden, output$params,c())
# Remove latent dimensions
th <- 0.03 #threshold to filter out latent features that are not significant
feat_toRemove <- which(sum(output$hidden$Z) < N*th) # filter features with insufficient number of obs. assigned
hidden <- remove_dims(output$hidden, feat_toRemove)
sorted_patterns<- get_feature_patterns_sorted(hidden$Z)
B_aux<-matrix(unlist(hidden$B),nrow=dim(hidden$B)[1],ncol=dim(Zp)[2],byrow=TRUE)
Kest <-dim(B_aux)[2]
Zp <- diag(Kest)
Zp[,1] <- 1 # bias active
Zp <- Zp[1:(min(5,Kest)),]
leg <-c('F0','F1', 'F2', 'F3', 'F4')
colours<-c()
GLFM_plotPatterns(data, output$hidden,output$params, Zp, list(leg,colours))


