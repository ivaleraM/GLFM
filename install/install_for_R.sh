Rscript -e "setwd('../src/Ccode')"
Rscript -e 'install.packages("Rcpp")'
Rscript -e 'require("Rcpp")'
Rscript -e 'install.packages("GLFM/src/Ccode/wrapper_R", repos = NULL, type="source")'
Rscript -e 'setwd("../src/GLFMR")'

