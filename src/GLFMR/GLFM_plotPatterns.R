GLFM_plotPatterns<-function(data,hidden,params,Zp,varargin){
  
D<-dim(hidden$B[[1]])[1]
#print(D)
#for(dd in 1:length(p$Results$idxD)){
  for(d in 1:D){
 # d<-p$Results$idxD[dd]
  if(data$C[d]=='g'|| data$C[d] == 'p' || data$C[d] == 'n'){
    idxs_nans <- which(is.nan(data$X[,d]))
    if(length(idxs_nans) > 0){
      idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
      mm <- min(data$X[idxs_nonnans,d])
      MM <- max(data$X[idxs_nonnans,d])
      h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
      h$density <- h$counts/sum(h$counts)
      plot(h,freq=FALSE, main = "histogram of x (proportions)")
      par(new=T)
    }
    mm <- min(data$X[,d])
    MM <- max(data$X[,d]) 
    h <- hist(data$X[,d], breaks=(mm-1):(MM+0.5))
    h$density <- h$counts/sum(h$counts)
    plot(h,freq=FALSE, main = "histogram of x (proportions)")
    par(new=T)
    print("Press return to continue")
  }
    #duda aqui
    pdf_val<-GLFM_computePDF(data,Zp,hidden,params,d)
    if(data$C[d]=='c'||data$C[d] == 'o' ){
      leges <- computeLeg(Zp,c())
      condition <- rep(as.character(leges) , ncol(pdf_val$pdf))
      sa <- stack(as.data.frame(t(pdf_val$pdf)))
      sa$ind<-condition
      sa$x <- rep(seq_len(ncol(pdf_val$pdf)), nrow(pdf_val$pdf))
      ggplot(sa, aes(fill=ind, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
      print("Press return to continue")
    }
    else if(data$C[d] == 'n'){
      h <-plot(pdf_val$xd)
     
    }
    else{
      h <-plot(pdf_val$xd)
      print("Nominal data")
    }
  }
}
  #To be completed
# Table 2 for prostate cancer types of data
idxs_nans <- which(is.nan(data_prost$X[,5]))
idxs_nonnans<-setdiff(1:(length(data_prost$X[,5])),idxs_nans)
   mm <- min(data_prost$X[idxs_nonnans,5])
   MM <- max(data_prost$X[idxs_nonnans,5])
   h <- hist(data_prost$X[idxs_nonnans,5], breaks=(mm-1):(MM+0.5))
   h$density <- h$counts/sum(h$counts)
   plotcols<-c('red','blue','green','pink','yellow')
   plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="b")
   par(new=T)
   for(d in 2:5){
plot(pdf_val$xd,pdf_val$pdf[d,],xlab = "",ylab="",col=plotcols[d],type="b",axes=F)
     par(new=T)
#print("Press return to continue")
   }
   #par(new=F)
plot(h,freq=FALSE, main = "histogram of x (proportions)",axes=F)
par(new=F)

# if(length(idxs_nans) > 0){
#   idxs_nonnans<-setdiff(1:(length(data_prost$X[,4])),idxs_nans)
#   mm <- min(data_prost$X[idxs_nonnans,d])
#   MM <- max(data_prost$X[idxs_nonnans,d])
#   h <- hist(data_prost$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
#   h$density <- h$counts/sum(h$counts)
#   plot(h,freq=FALSE, main = "histogram of x (proportions)")
# }