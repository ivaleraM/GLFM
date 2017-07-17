
f_p_1<-function(x,mu,w){
  if (w == 0){
    stop( 'scaling factor should never be 0')
  }
  else{
    y <- log( exp(w*(x-mu) - 1) )  
  }
}



