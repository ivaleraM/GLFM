#include <RcppGSL.h> 
// [[Rcpp::depends(RcppGSL)]]
//#include <gsl/gsl_sf_exp.h>
//#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_randist.h>
using namespace Rcpp;
// [[Rcpp::export]]
#include "../core/GeneralFunctions.cpp"
#include "../core/InferenceFunctions.cpp"
#include <InferenceFunctions.h> 
RcppExport SEXP initialize_wrapper(SEXP N_, SEXP D_, SEXP maxK_, SEXP missing_, SEXP X_, SEXP C_, SEXP B_, SEXP theta_, SEXP R_, SEXP f_, SEXP mu_,  SEXP w_, SEXP s2Y_){
 
 // Convert the input to C++ types
 int N = as<int>(N_), D = as<int>(D_), maxK = as<int>(maxK_), R = as<int>(R_);
double missing = as<double>(missing_), f = as<double>(f_), mu = as<double>(mu_), w = as<double>(w_),  s2y = as<double>(s2y_);
char C = as<char>(C_);
RcppGSL::matrix<double> X = gsl_matrix(X_), B =  gsl_matrix(B_); 
RcppGSL::vector<double> theta = gsl_vector(theta_);
//RcppGSL::vector_view<double> colview = gsl_matrix_const_column(G, j);

// call the underlying C++ function
int maxR = initialize_func (N, D, maxK, missing,X, C, B, theta, R, f, mu, w, s2Y)
  
// return the result as SEXP
return wrap(maxR);
}

RCPP_MODULE(mod){
  function("maxR", &maxR, "Provides some number");
}  





