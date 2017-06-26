#' @description: transformation function for real-valued data
#' Y -> X (from pseudo-obversations to data)
#' @param y: N*R where N is the number of observations and R is the num of categories
#' @param mu: mean to shift data
#' @param w: weights
#' @return shifted data

f_g<-function(y,mu,w){
  if(w==0){
    stop('scaling factor should never be 0') 
  }
  else{
    x<-y/w+mu
  }
}


