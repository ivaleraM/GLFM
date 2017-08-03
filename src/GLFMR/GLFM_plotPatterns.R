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
      sa$Patterns <-condition
      specie <-plotlabels[[d]]
      sa$x <-specie
     # x11()
      p1<-ggplot(sa, aes(fill=Patterns, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
      plot.new()
      print(p1)
      #print(sa)
      print("Press return to continue")
    }
    else if(data$C[d] == 'n'||data$C[d] == 'p'||data$C[d] == 'g'){
      # As many colours as the number of patterns to plot P
      plot.new()
      idxs_nans <- which(is.nan(data$X[,d]))
      idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
      mm <- min(data$X[idxs_nonnans,d])
      MM <- max(data$X[idxs_nonnans,d])
      if(params$transf_dummie){
        if(params$transf_dummie && d == params$idx_transform){
          print("entra aqui")
         # mm <- params$t_1(mm)
        #  MM <- params$t_1(MM)
          h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
          h$density <- h$counts/sum(h$counts)
          plot(h,freq=FALSE, main = "histogram of x (proportions)",axes=F)
          par(new=T)
          plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="l")
          par(new=T)
          for(pp in 2:P){
            plot(pdf_val$xd,pdf_val$pdf[pp,],xlab = "",ylab="",col=plotcols[pp],type="l",axes=F)
            par(new=T)
            #  # print("Press return to continue")
          }
          par(new=F)
        }
      }else{
       plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="l")
       par(new=T)
       for(pp in 2:P){
         plot(pdf_val$xd,pdf_val$pdf[pp,],xlab = "",ylab="",col=plotcols[pp],type="l",axes=F)
         par(new=T)
      #  # print("Press return to continue")
       }
       par(new=F)
      }
    }
  }
}

# Table 2 for prostate cancer types of data

   
