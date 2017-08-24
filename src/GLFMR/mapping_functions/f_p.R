#' @description Transformation function for positive data
#' @param y: N*K-sized matrix, where N is the number of observations, K categories
#' @param mu: mean to shift data
#' @param w: weights
#' @Return x transformed data


f_p<-function(y,mu,w){
 if(w == 0){
    stop('scaling factor should never be 0')
 }
  else{
    x <- log( exp(y) + 1 )/w + mu
  }
  return(x)
}