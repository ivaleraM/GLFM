#var1<-c("HOLA")
# sapply(var1, tolower)
#varargin is a list
#' @return A list with posterior draws

GLFM_infer<-function(data,varargin){
  varargin_size<-length(varargin)
  if(varargin_size==0){
    hidden<-c()
    params <- init_default_params(data, c())
  }
  else if(varargin_size<3)
  {
    hidden<-varargin[1]
    switch(varargin_size,params <- init_default_params(data, c()) ,params <- init_default_params(data, unlist(varargin[2])))
  }
  else
  {
    stop("Incorrect number of input parameters: should be 0, 1 or 2")
  }
  D <- dim(data$X)[2]
  N <- dim(data$X)[1]
  if(length(hidden)==0)
  {
    m0<-matrix(0,N,2)
    Z <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.8,0.2)))
  }
  if(params$bias == 1){
    Z <-cbind(rep(1,N),Z)
  }
  else if(params$bias>1){
    stop("There is more than 1 bias specified, but the structure Z has not been initialized") 
  }
  # replace missing values
  data$X[which(is.nan(data$X))] <- params$missing
  # Change labels of categorical or ordinal data
  # Cambiamos las etiquetas de los datos categoricos o ordinales que no son datos
  # faltantes para que su respectiva categoria empiece en uno. y por eso se toma
  # el minimo, etc.
  # ---To be completed ---
  # From the posterior
 R<-rep(1,D)
  return(list("Z"=Z,"B"=B,"theta"=theta,"mu"= mu,"w"=w,"s2Y"=s2y,"R"=R))
}
 
  
  
  
  