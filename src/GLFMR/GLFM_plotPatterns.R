GLFM_plotPatterns<-function(data,hidden,params,Zp,varargin){
leges<-varargin$leges
plotcols<-varargin$colours
D<-dim(hidden$B)[1]
P<-dim(Zp)[1]
  for(d in 1:D){
    x11()
    pdf_val<-GLFM_computePDF(data,Zp,hidden,params,d)
    idxs_nans <- which(is.nan(data$X[,d]))
    idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
    mm <- min(data$X[idxs_nonnans,d])
    MM <- max(data$X[idxs_nonnans,d])
    h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
    if(data$C[d]=='c'| data$C[d] == 'o' ){
      # Adds the empirical
      auxpdf<-rbind(h$counts/sum(h$counts),pdf_val$pdf)
      condition <- rep(as.character(leges) , ncol(auxpdf))
      sa <- stack(as.data.frame((auxpdf)))
      sa$Patterns <-condition
      specie <-plotlabels[[d]]
      sa$x <-specie
      p1<-ggplot(sa, aes(fill=Patterns, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
      p1<-p1+ ggtitle(plottitles[[d]])
      plot.new()
      print(p1)
      print("Press return to continue")
    }
    else if(data$C[d] == 'n'||data$C[d] == 'p'||data$C[d] == 'g'){
      plot.new()
        if(params$transf_dummie && d == params$idx_transform){
          h$density <- h$counts/sum(h$counts)
          plot(h,xlab=expression("x"[d]), ylab=expression("pdf"[x]),freq=FALSE, main = plottitles[[d]],axes=F)
          par(new=T)
          plot(pdf_val$xd,pdf_val$pdf[1,],xlab = expression("x"[d]),ylab=expression("pdf"[x]),col=plotcols[1],type="l")
          par(new=T)
          for(pp in 2:P){
            plot(pdf_val$xd,pdf_val$pdf[pp,],xlab = expression("x"[d]),ylab=expression("pdf"[x]),col=plotcols[pp],type="l",axes=F)
            par(new=T)
          }
          par(new=F)
        }
      else{
        h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
        h$density <- h$counts/sum(h$counts)
        plot(h,freq=FALSE,xlab=expression("x"[d]), ylab=expression("pdf"[x]), main = plottitles[[d]],axes=F)
        par(new=T)
       plot(pdf_val$xd,pdf_val$pdf[1,],xlab = expression("x"[d]), ylab=expression("pdf"[x]),col=plotcols[1],type="l")
       par(new=T)
       for(pp in 2:P){
         plot(pdf_val$xd,pdf_val$pdf[pp,],xlab =expression("x"[d]), ylab=expression("pdf"[x]),col=plotcols[pp],type="l",axes=F)
         par(new=T)
       }
       par(new=F)
      }
    }
  }
}


   
