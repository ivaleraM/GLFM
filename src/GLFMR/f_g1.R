
#' @description  Transformation function for real-valued data
#' Y -> X (from pseudo-obversations to data)
#' @param y: N*R where N is the number of observations and R is the num of categories
#' @param mu: mean to standarise data
#' @param w: weights
#' @Return x standarised data

f_g1<function(x, mu, w){
  if (w == 0){
    stop('scaling factor should never be 0')
  }
  else{
    y <- w * (x - mu)
  }
}

