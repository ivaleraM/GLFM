#' @param y: is a column vector
#' @param Thet: is a column vector too that is used for binning
#' 
f_o<-function(y,Thet){
  x<-rep(0,length(y))
  for(j in 1:length(Thet)){
    if(j == 1){
      mask<-which(y<=Thet[1])
    }
    else{
      mask<-which(y>Thet[j-1]&& y<=Thet[j] )
    }
    x[mask]<-j
  }
  x [which( x == 0)] <- length(Thet)+1
  return(x)
}