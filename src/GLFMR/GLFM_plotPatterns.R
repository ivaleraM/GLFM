GLFM_plotPatterns<-function(data,hidden,params,Zp,varargin){
leges<-varargin$leges
plotcols<-varargin$colours
D<-dim(hidden$B)[1]
P<-dim(Zp)[1]
#print(D)
#for(dd in 1:length(p$Results$idxD)){
  for(d in 1:D){
    x11()
    pdf_val<-GLFM_computePDF(data,Zp,hidden,params,d)
    if(data$C[d]=='c'| data$C[d] == 'o' ){
     # leges <- computeLeg(Zp,c())
      condition <- rep(as.character(leges) , ncol(pdf_val$pdf))
      sa <- stack(as.data.frame((pdf_val$pdf)))
      sa$ind<-condition
      sa$x <- rep(seq_len(ncol(pdf_val$pdf)), nrow(pdf_val$pdf))
     # x11()
      p1<-ggplot(sa, aes(fill=ind, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
      plot.new()
      print(p1)
      print(sa)
      print("Press return to continue")
    }
    else if(data$C[d] == 'n'||data$C[d] == 'p'){
      # As many colours as the number of patterns to plot P
      plot.new()
      idxs_nans <- which(is.nan(data$X[,d]))
      idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
      mm <- min(data$X[idxs_nonnans,d])
      MM <- max(data$X[idxs_nonnans,d])
      h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
      h$density <- h$counts/sum(h$counts)
      plot(h,freq=FALSE, main = "histogram of x (proportions)",axes=F)
      par(new=T)
       plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="b")
       par(new=T)
       for(d in 2:P){
         plot(pdf_val$xd,pdf_val$pdf[d,],xlab = "",ylab="",col=plotcols[d],type="b",axes=F)
         par(new=T)
      #  # print("Press return to continue")
       }
       par(new=F)
      
      #duda aqui
      #par(new=F)
      #Includes the histogram from above
      #plot(h,freq=FALSE, main = "histogram of x (proportions)")
      #par(new=T)
      
    }
    else{
     # h <-plot(pdf_val$xd)
      print("positive or real data")
    }
    # if(data$C[d]=='g'|| data$C[d] == 'p' || data$C[d] == 'n'){
    #   idxs_nans <- which(is.nan(data$X[,d]))
    #   idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
    #   mm <- min(data$X[idxs_nonnans,d])
    #   MM <- max(data$X[idxs_nonnans,d])
    #   h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
    #   h$density <- h$counts/sum(h$counts)
    #   plot(h,freq=FALSE, main = "histogram of x (proportions)",axes=F)
    #  par(new=F)
    # }
  }
}

# Table 2 for prostate cancer types of data

   
