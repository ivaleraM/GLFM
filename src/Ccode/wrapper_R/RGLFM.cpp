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
#include <gsl/gsl_randist.h>
using namespace Rcpp;
// [[Rcpp::export]]
//#include "../core/GeneralFunctions.cpp"
#include "../core/InferenceFunctions.cpp"
#include "../core/InferenceFunctions.h"
#include ".../core/GeneralFunctions.h"
// #include <InferenceFunctions.h> 
RcppExport SEXP initialize_wrapper(SEXP N_, SEXP D_, SEXP K_, SEXP maxK_, SEXP missing_, const RcppGSL::Matrix & X, const RcppGSL::Matrix & Z, SEXP C_, const RcppGSL::Matrix & B, const RcppGSL::Vector & theta, SEXP R_, SEXP f_, SEXP mu_,  SEXP w_, SEXP s2y_, SEXP bias_, SEXP alpha_, SEXP s2B_, SEXP s2u_, SEXP Nsim_){
 
 // Convert the input to C++ types
 int N = as<int>(N_), D = as<int>(D_), maxK = as<int>(maxK_), R = as<int>(R_), Nsim = as<int>(Nsim_), K = as<int>(K_);
double missing = as<double>(missing_), f = as<double>(f_), mu = as<double>(mu_), w = as<double>(w_),  s2Y = as<double>(s2y_); 
  double bias = as<double>(bias_), alpha = as<double>(alpha_), s2B = as<double>(s2B_), s2u = as<double>(s2u_) ; 
char C = as<char>(C_);

// call the underlying C++ function
int maxR[] = initialize_func (N, D, maxK, missing,X, C, B, theta, R, f, mu, w, s2Y);
int bu[] = IBPsampler_func (missing, X, C, Z, B, theta, R, f, mu, w, maxR, bias, N, D, K, alpha, s2B, s2Y,  s2u, maxK, Nsim);
  
  print "maxR = %d" % maxR;
  //{
// return the result as SEXP
//return wrap(maxR);
}
//RCPP_MODULE(mod){
//  function("maxR", &maxR, "Provides some number");
//}  