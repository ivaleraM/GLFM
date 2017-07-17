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
sa <- stack(as.data.frame(t(pdf_val$pdf)))
sa$ind<-condition
sa$x <- rep(seq_len(ncol(pdf_val$pdf)), nrow(pdf_val$pdf))
ggplot(sa, aes(fill=ind, y=values, x=x)) + geom_bar(position = "dodge", stat="identity") 
#+    facet_wrap(~ind)


