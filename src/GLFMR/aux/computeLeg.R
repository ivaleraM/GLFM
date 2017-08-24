# num2str(x) en Matlab es as.character(x)

computeLeg<-function(Zp,params){
  if(length(params$bias)>0){
    Zp <-Zp[,-1]
    n_leg<-dim(Zp[,-1])[1]
  }
  else{
   n_leg<- dim(Zp)[1]
  }
  legends<-list()
  for(jj in 1:n_leg){
    print(paste(as.character(Zp[jj,]),collapse =" ")) 
    legends<-c(legends, paste(as.character(Zp[jj,]),collapse =" "))
  }
  return(legends)
}
