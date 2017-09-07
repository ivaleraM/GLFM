#' @description A function that plots the pdf/probability functions
#' Inputs:
#'@param data is a list with X and C
#'@param X: N*D data matrix X
#'@param C: 1xD string with data types, D = number of dimensions
#'@param Zp: P*K matrix of patterns (P is the number of patterns)
#'@param hidden: list with the following latent variables learned by the model (B,Z):
#'@param B: latent feature list with D matrices of size  K * maxR  where
#'@param D: number of dimensions
#'@param K: number of latent variables
#'@param maxR: maximum number of categories across all dimensions
#'@param Z: PxK matrix of feature patterns
#'@param params: list with parameters (mu,w,theta)
#'@param  mu: 1*D shift parameter
#'@param w:  1*D scale parameter
#'@param  theta: D*maxR matrix of auxiliary vars (for ordinal variables)
# ----------------(optional) ------------------
#'@param varargin a list with extra arguments: list of legends and colours (leges,plotcols)
#'@param leges string of legends for patterns
#'@param plotcols vector of colours

GLFM_plotPatterns<-function(data,hidden,params,Zp,varargin){
  if(length(varargin)>0){
#leges<-varargin$leges
leges<-varargin[[1]]
#plotcols<-varargin$colours
plotcols<-varargin[[2]]
  }
D<-dim(hidden$B)[1]
P<-dim(Zp)[1]
  for(d in 1:D){
    x11()
    pdf_val<-GLFM_computePDF(data,Zp,hidden,params,d)
    idxs_nans <- which(is.nan(data$X[,d]))
    idxs_nonnans<-setdiff(1:(length(data$X[,d])),idxs_nans)
    mm <- min(data$X[idxs_nonnans,d])
    MM <- max(data$X[idxs_nonnans,d])
    if(data$C[d]=='c'| data$C[d] == 'o' ){
      # Adds the empirical
      h <- hist(data$X[idxs_nonnans,d], breaks=(mm-1):(MM+0.5))
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
    }
    else if(data$C[d] == 'n'||data$C[d] == 'p'||data$C[d] == 'g'){
      plot.new()
          if(params$transf_dummie && d %in% params$idx_transform){
          h <- hist(data$X[idxs_nonnans,d], breaks=100,prob=TRUE,xlab=expression("x"[d]), ylab=expression("pdf"[x]),main = plottitles[[d]])
          lines(pdf_val$xd,pdf_val$pdf[1,],xlab = expression("x"[d]),ylab=expression("pdf"[x]),col=plotcols[1],type="l")

          for(pp in 2:P){
            lines(pdf_val$xd,pdf_val$pdf[pp,],xlab = expression("x"[d]),ylab=expression("pdf"[x]),col=plotcols[pp],type="l")
          }
        }
      else{
        plot.new()
        h <- hist(data$X[idxs_nonnans,d], breaks=100,prob=TRUE,xlab=expression("x"[d]), ylab=expression("pdf"[x]), main = plottitles[[d]])
       lines(pdf_val$xd,pdf_val$pdf[1,],xlab = expression("x"[d]), ylab=expression("pdf"[x]),col=plotcols[1],type="l")
       for(pp in 2:P){
         lines(pdf_val$xd,pdf_val$pdf[pp,],xlab =expression("x"[d]), ylab=expression("pdf"[x]),col=plotcols[pp],type="l")
       }
      }
    }
  }
}


   
