#' Function to initialize or complete params structure with all
#' simulation parameters and hyperparameters for the GLFM
#' @param  data is a list with X, Z, C in order to be able to initialize
#' @param params with data-based values
#' @param t: eventual external transform of obs. X = params.t{d}(Xraw)
#' @param  t_1: inverse transform
#' @param dt_1:derivative of the inverse transform

init_default_params<-function(data,params)
{
params_names<-c("missing","alpha","bias","s2u","s2B","Niter","maxK","verbose","numS","t","t_1","dt_1")
idx_to_fill<-which(params_names %in% params)

filled_param_names<-paste("params",params_names[idx_to_fill],sep=".")
param_values<-list(-1,1,0,0.01,1,1000,dim(data)[2],1,100,list(1,dim(data$X)[2]),list(1,dim(data$X)[2]),list(1,dim(data$X)[2]))
#names(param_values)<-param_names
params_to_return<-param_values[idx_to_fill]
names(params_to_return)<-filled_param_names
return(params_to_return)
}
