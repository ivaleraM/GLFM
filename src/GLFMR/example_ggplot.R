# Pretty plots using ggplot2 for the barplots
require(ggplot2)
# First run demo_data_exploration_prostate.R until before the plotpatterns routine
# dimensions 1 and 3 are categorical
pdf_val<-GLFM_computePDF(data_prost,Zp,output$hidden,output$params,3)
Kest <-dim(output$hidden$B[[1]])[1]
Zp <- diag(Kest)
Zp[,1] <- 1 # bias active
Zp <- Zp[1:(min(5,Kest)),]
leges<-computeLeg(Zp,c())
condition=rep(as.character(leges) , ncol(pdf_val$pdf))
sa <- stack(as.data.frame((pdf_val$pdf)))
sa$ind<-condition
sa$x <- rep(seq_len(ncol(pdf_val$pdf)), nrow(pdf_val$pdf))
ggplot(sa, aes(fill=ind, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
#+    facet_wrap(~ind)

# problema con los datos ordinales, revisar
 pdf_val<-GLFM_computePDF(data_prost,Zp,output$hidden,output$params,2)
# 
# leges <- computeLeg(Zp,c())
# condition <- rep(as.character(leges) , ncol(pdf_val$pdf))
 sa <- stack(as.data.frame((pdf_val$pdf)))
 sa$ind<-condition
 sa$x <- rep(seq_len(ncol(pdf_val$pdf)), nrow(pdf_val$pdf))
 ggplot(sa, aes(fill=ind, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 

# Mejorar la grafica de histograma y probability mass function
 idxs_nans<-which(data_prost$X[,4])
 idxs_nonnans<-setdiff(1:(length(data_prost$X[,4])),idxs_nans)
 mm <- min(data_prost$X[idxs_nonnans,4])
 MM <- max(data_prost$X[idxs_nonnans,4])
 h <- hist(data_prost$X[idxs_nonnans,4], breaks=(mm-1):(MM+0.5))
 h$density <- h$counts/sum(h$counts)
 pdf_val<-GLFM_computePDF(data_prost,Zp,output$hidden,output$params,4)
 plotcols<-c('red','blue','green','pink','yellow')
 plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="b")
 par(new=T)
 for(d in 2:5){
   plot(pdf_val$xd,pdf_val$pdf[d,],xlab = "",ylab="",col=plotcols[d],type="b",axes=F)
   par(new=T)
   #print("Press return to continue")
 }
 #par(new=F)
 #Includes the histogram from above
 plot(h,freq=FALSE, main = "histogram of x (proportions)",axes=F)
 par(new=F)
 plot(h,freq=FALSE, main = "histogram of x (proportions)")
 par(new=T)
 plot(pdf_val$xd,pdf_val$pdf[1,],xlab = "",ylab="",col=plotcols[1],type="b")
 

pdf_val<-GLFM_computePDF(data_prost,Zp,output$hidden,output$params,5)
idxs_nans <- which(is.nan(data_prost$X[,5]))
 idxs_nonnans<-setdiff(1:(length(data_prost$X[,5])),idxs_nans)
 mm <- min(data_prost$X[idxs_nonnans,5])
 MM <- max(data_prost$X[idxs_nonnans,5])
 h <- hist(data_prost$X[idxs_nonnans,5], breaks=(mm-1):(MM+0.5))
 h$density <- h$counts/sum(h$counts)
 leges <- computeLeg(Zp,c())
 condition <- rep(as.character(leges) , ncol(pdf_val$pdf))
 sa_den <- stack(as.data.frame(t(pdf_val$pdf)))
 #sa_den$ind<-condition
 sa_den$x<-rep(pdf_val$xd,nrow(pdf_val$pdf))
# Density plots per patterns
 p9<-ggplot(sa_den, aes(x=x)) +  geom_line(aes(y=values, color=ind))
 p9 + theme(legend.position="none")
 p9 + scale_fill_discrete(labels=c(unlist(leges)))
 #adding the right legend
# ugly histogram here
# dataset2<-data.frame("yy"=h$counts,"xx"=h$density)
# p8 <-ggplot(dataset2,aes(x=h$counts))+geom_histogram(bins=100,binwidth = 5)