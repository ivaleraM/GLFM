#' transformation function for categorical data
#' @param input : is the label
#' @param y: N*K-sized matrix, where N is the number of observations, K categories
#'@import matrixStats for colMaxs
#'@return a: the values of the maximum element per row and the index of the first max in x
f_c<-function(y){
  a<-max.col(y)
  #if(length(x)>1){
  #  x<-X[1]
  #}
  return(list(a))
# [a,x] = max(y,[],2)
}


