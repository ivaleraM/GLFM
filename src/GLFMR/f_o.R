#' @param y: N*K-sized matrix, where N is the number of observations, K categories
#' @param Thet: is a column vector too that is used for binning
#' 
f_o<-function(y,Thet){
  x<-rep(0,length(y))
  for(j in 1:length(Thet)){
    if(j == 1){
      mask<-(y<=Thet[1])
    }
    else{
      mask<-(y>Thet[j-1]&& y<=Thet[j] )
    }
    x[mask]<-j
  }
  x [( x == 0)] <- length(Thet)+1
  return(x)
}