# demo_data_exploration_counties.R
rm(list=ls())
graphics.off()
require(R.matlab)
require(ggplot2)
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("GLFM_computePDF.R")
source("GLFM_plotPatterns.R")
#source("remove_dims.R")
datos_counties<-readMat('counties.mat')
source("aux/get_feature_patterns_sorted.R")
source("aux/computeLeg.R")
Xauxi <- as.matrix(datos_counties$data[1,1,1][[1]],ncol=19,nrow= 3141, byrow=TRUE)
Xfull<-matrix(Xauxi,nrow=3141,ncol=19)
Xfull[,9]<-Xfull[,8]+Xfull[,9]
C<-unlist(datos_counties$data[2,1,1],use.names = FALSE)
Cfull<-strsplit(as.character(C), "")
cat_labels_full <-unlist(datos_counties$data[3,1,1],use.names = FALSE)
y_labels_full<-unlist(datos_counties$data[4,1,1],use.names = FALSE)
idx_to_remove <- c(1,3,4, 6, 7,8, 10,19) 
idx_to_keep <- setdiff(1:19,idx_to_remove)
X<-Xfull[,idx_to_keep]
C<-Cfull[[1]][idx_to_keep]
N<-dim(X)[1]
D<-dim(X)[2]
param_names<-c("missing","s2u","s2B","alpha","Niter","maxK","bias","transf_dummie")
missing <- -1
s2Y <- 0   
s2u <- .005 
s2B <- 1   
alpha <- 1   
Niter <- 100
maxK <- 10;
bias <- 1
func <- 1*rep(1,D)
transf_dummie <-TRUE
if(transf_dummie){
  # this should be a list of lenght equal to the number of different transformations
  idx_transform1 <- c(2, 5, 6, 11)
  idx_transform2 <- 10
  idx_transform <- list(idx_transform1,idx_transform2)
  # transformation to apply to raw data
  t_11<-function(x){log(x+1)}
  # inverse transform to recover raw data
  t_inv1<-function(y){exp(y)-1}
  # derivative of transform
  dt_11<-function(x){1/as.vector(x+1)}
  t_12<-function(x){log((100-x)+1)}
  # inverse transform to recover raw data
  t_inv2<-function(y){-exp(y)+101}
  # derivative of transform
  dt_12<-function(x){-1/as.vector(101-x)}
  # change type of data due to transformation
  ext_datatype <-list('p','p')
  t_1<-list(t_11,t_12)
  dt_1<-list(dt_11,dt_12)
  t_inv<-list(t_inv1,t_inv2)
  param_names<-c(param_names,'t_1','dt_1','t_inv','ext_datatype','idx_transform')
  params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias,transf_dummie,t_1,dt_1,t_inv,ext_datatype,idx_transform)
 }else{ 
  params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias,transf_dummie)
    }
names(params)<-param_names
m0<-matrix(0,nrow=N,ncol=1)
Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.2,0.8)))
data_counties<-list("X"=X,"C"=C)
output <- GLFM_infer(data_counties, list(Z,params))
X_map <- GLFM_computeMAP(data$C, output$hidden$Z, output$hidden, output$params,c())
# To be completed

