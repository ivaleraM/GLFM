# Tiny example of how to use the Rcpp package to call C++ files
# If you don't have the packages do install.packages("Rcpp")
require(Rcpp)
require(RcppEigen)
require(RcppGSL)# For all GSL library bits
setwd("~/Documents/Working_papers/FAP_Rpackage/GLFM/Rpackage_create/Rcpp_examples/")
# sourceCpp(".cpp") is the equivalent of source for C++
sourceCpp("mvnorm.cpp")
# Generates a sample from a multivariate Normal distribution
# Note that the name of the function is different than the name of the file 
# which is not great. I always name my files as my functions and include just
# one function per file so the code is modular.
rmvnorm(10,matrix(10,nrow=2,ncol=1),matrix(c(c(1,0),c(0,1)),nrow=2,ncol=2))
