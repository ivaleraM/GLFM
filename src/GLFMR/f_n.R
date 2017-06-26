#' @description  Transformation function for count data
#' @param y: N-sized vector, where N is the number of observations 
#' @param mu: mean to shift data
#' @param w: weights
#' @Return x transformed data


f_n<-function(y,mu,w){
  if (w == 0){
    stop('scaling factor should never be 0')
  }
  else{
    x <- floor(log( exp(y) + 1 )/w + mu )
  }
    return(x)
}


