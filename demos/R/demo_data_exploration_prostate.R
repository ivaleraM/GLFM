# demo_data_exploration_prostate
rm(list=ls())
graphics.off()
require(R.matlab)
require(ggplot2)
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("GLFM_computePDF.R")
source("GLFM_plotPatterns.R")
source("aux/remove_dims.R")
datos_prostate<-readMat('prostate_v3.mat')
source("aux/get_feature_patterns_sorted.R")
source("aux/computeLeg.R")
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
#idx_cat1<-which(X[,1]==3)
#idx_cat2<-which(X[,1]==4)
#Xauxcat<-X[,1]
#Xauxcat[idx_cat1]<-1
#Xauxcat[idx_cat2]<-2
#aux11<-rep("stage",dim(X)[1])
#cat_labels1 <-paste(aux11,Xauxcat)
aux1<-rep(paste("stage",cat_labels_full[1:2]),6)
cat_labels1<-aux1[order(aux1)]
aux2<-rep(paste(cat_labels_full[3:5], "mg"),6)
ord_labels2<-aux2[order(aux2,decreasing=FALSE)]
aux3<-rep(cat_labels_full[6:9],6)
cat_labels3<-aux3[order(aux3,decreasing=FALSE)]
y_labels<-y_labels_full[idx_toKeep]
y_labels_long<-y_labels_long_full[idx_toKeep]
N<-dim(X)[1]
D<-dim(X)[2]
plotlabels<-list(cat_labels1,ord_labels2,cat_labels3)
plottitles<-list("Type of cancer","Prognosis status", "Drug level", "Size of primary tumor (cm^2)","Serum prostatic acid phosphatase")
# pre-transform a subset of variables
#params
param_names<-c("missing","s2u","s2B","alpha","Niter","maxK","bias","transf_dummie","plotlabels","plottitles")
missing<--1
s2u<-0.005
s2B<-1
alpha<-1
Niter<-1000
maxK<-10
bias<-1
transf_dummie <-TRUE
  if(transf_dummie){
    idx_transform <- D # we transform the last dimension
    # transformation to apply to raw data
    t_1<-function(x){log(x+1)}
    # inverse transform to recover raw data
    t_inv<-function(y){exp(y)-1}
    # derivative of transform
    dt_1<-function(x){1/as.vector(x+1)}
    # change type of data due to transformation
    ext_datatype <-'p'
  param_names<-c(param_names,'t_1','dt_1','t_inv','ext_datatype','idx_transform')
  params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias,transf_dummie,plotlabels,plottitles,t_1,dt_1,t_inv,ext_datatype,idx_transform)
    } else{ 
      params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias,transf_dummie,plotlabels,plottitles)
          }
names(params)<-param_names

# Inference
Z<-c()
data_prost<-list("X"=X,"C"=C)
output <- GLFM_infer(data_prost, list(Z,params))
#Predict MAP estimate for the whole matrix X
X_map <- GLFM_computeMAP(data_prost$C, output$hidden$Z, output$hidden, output$params,c())
# Remove latent dimensions
th <- 0.03 #threshold to filter out latent features that are not significant
feat_toRemove <- which(sum(output$hidden$Z) < N*th) # filter features with insufficient number of obs. assigned
if(length(feat_toRemove)>0){
  output$hidden <- remove_dims(output$hidden, feat_toRemove)
  sorted_patterns <- get_feature_patterns_sorted(hidden$Z,c())
}
sorted_patterns<- get_feature_patterns_sorted(output$hidden$Z,c())
Kest <-dim(output$hidden$B[[1]])[1]
Zp <- diag(Kest)
Zp[,1] <- 1 # bias active
Zp <- Zp[1:(min(5,Kest)),]
leges <- computeLeg(rbind(rep(0, ncol(Zp)),Zp),c())
colours<-c('red','blue','green','pink','yellow')
# Falta calcular la probabilidad empirica (es lo que llaman baseline)

GLFM_plotPatterns(data_prost,output$hidden,output$params,Zp, list("leges"=leges,"colours"=colours) )


