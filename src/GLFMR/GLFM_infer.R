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
  # From the posterior
  return(list("Z"=hidden))
}
 
  
  
  
  