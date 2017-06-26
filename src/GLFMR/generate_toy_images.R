
#' @param s2x: observation noise
#' @param N: number of samples
#' @return data: data structure with X: N*D obs matrix and
#' @return  C: 1*D datatype string vector
%

 generate_toy_images<-function(N,s2x){
Btrue<-matrix(c(c(0,1.0,0,0,0,0,1,1,1,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),
                c(0,0.0,0,1,1,1,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),
              c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0, 1,1,1,0,0,0),
              c(0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,1,1,1, 0,0,0,0,1,0, 0,0,0,0,1,0)),
              nrow=4,ncol=36,byrow=TRUE)
K <- dim(Btrue)[1] # number of latent features
D <- dim(Btrue)[2] # number of dimensions
m0<-matrix(0,nrow=N,ncol=K)
Ztrue <- apply(m0, c(1,2), function(x) sample(c(0,1),1,prob=c(0.2,0.8)))
norm_mat<-matrix(rnorm(N*D),nrow=N,ncol=D,byrow=TRUE)
X <- sqrt(s2x) * norm_mat +  Ztrue %*% Btrue 
C <- rep('g',D)
data<-list("X"=X,"C"=C)
gT<-list("B"=Btrue,"Z"=Ztrue)
return(list(data,gT))
}

