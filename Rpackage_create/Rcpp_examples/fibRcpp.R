# First example using inline

incltxt<-'
int fibonacci(const int x){
if (x == 0) return(0);
if (x == 1) return(1);
return fibonacci(x-1)+fibonacci(x-2);
}'

fibRcpp <- cxxfunction(signature(xs="int"), plugin="Rcpp",incl = incltxt, body = 'int x = Rcpp::as<int>(xs);
                       return Rcpp::wrap(fibonacci(x) );')
