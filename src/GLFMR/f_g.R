#' @description: transformation function for real-valued data
#' Y -> X (from pseudo-obversations to data)
#' @param y: N*K-sized matrix, where N is the number of observations, K categories
#' @param mu: mean to shift data
#' @param w: weight to scale it
#' @Return x transformed data

f_g<-function(y,mu,w){
  if(w==0){
    stop('scaling factor should never be 0') 
  }
  else{
    x<-y/w+mu
  }
}


