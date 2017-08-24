# demo_data_exploration_counties.R
rm(list=ls())
graphics.off()
require(R.matlab)
require(ggplot2)
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("GLFM_computePDF.R")
source("GLFM_plotPatterns.R")
source("remove_dims.R")
source("aux/get_feature_patterns_sorted.R")
source("aux/computeLeg.R")
datos_counties<-readMat('../../datasets/counties.mat')

# Missing dataset file and data processing here

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
  idx_transform <- D 
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
  ext_datatype <-'p'
  param_names<-c(param_names,'t_11','dt_11','t_inv1','t_12','dt_12','t_inv2','ext_datatype','idx_transform')
  params<-list(missing,s2u,s2B,alpha,Niter,maxK,bias,transf_dummie,t_11,dt_11,t_inv1,,t_12,dt_12,t_inv2,ext_datatype,idx_transform)
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

