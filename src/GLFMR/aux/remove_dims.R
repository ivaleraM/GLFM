#' Removes dimensions
remove_dims<-function(hidden,idxs){
  hidden$Z[,idxs]<-c()
  hidden$B[,idxs]<-c()
  return(hidden)
}