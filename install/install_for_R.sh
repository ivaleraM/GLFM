#!/bin/bash

Rscript -e 'install.packages("Rcpp")'
Rscript -e 'install.packages("RcppGSL")'
Rscript -e 'install.packages("matrixStats")'
Rscript -e 'install.packages("ggplot2")'
Rscript -e 'install.packages("R.matlab")'

cd ../src/Ccode/wrapper_R
R CMD build RcppGSLExample
R CMD INSTALL RcppGSLExample

