#' @description A function to return matrix with all patterns dim = (numPatterns*numFeatures)
#' @param  Z: N*K matrix
#' @param fid: file identifier where to print description (default: print in command window with fid=1)
#' @param patterns: P*K (complete list of patterns, this might be  useful if we already know the list of patterns and Z is a subset of the
#' matrix, so not all patterns might be present in that subset).
#' @return patterns: (numPatterns*numFeatures) matrix
#' @return C is an assignment vector, with the pattern id for each patient (numPatients*1)

get_feature_patterns_sorted<-function(Z,varargin){
  N <- dim(Z)[1]
  D <- dim(Z)[2] # number of features
  C_types <- rep(0,N)
  if(length(varargin)==0){
    fid<-1
    patterns<-sapply(Z,unique)
    patterns<-Z[!duplicated(Z),]
  }
  else if(length(varargin)==1){
    fid <- varargin[[1]]
    patterns<-sapply(Z,unique)
  }
  else{
    fid <- varargin[[1]]
    patterns <- varargin[[2]]
  }
  L <- rep(0,dim(patterns)[1])
  for(r in 1:dim(patterns)[1] ){
    idxs_eq<-which(apply(Z, 1, function(x) all.equal(x , patterns[r,])) == "TRUE")
    C_types[idxs_eq]<-r
    L[r]<-length(idxs_eq)
  }
  L<-L[order(L,decreasing=TRUE)]
  patterns<-patterns[order(L,decreasing=TRUE),]
  return(list("paterns"=patterns,"types"=C_types))
}