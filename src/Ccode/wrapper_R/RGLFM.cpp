// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h> 
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

using namespace Rcpp;
// [[Rcpp::export]]
#include "../core/GeneralFunctions.cpp"
#include "../core/InferenceFunctions.cpp

List GLFM(NumericMatrix X, CharacterVector C, NumericMatrix Z, NumericMatrix B, double theta, int R, NumericVector weights, int maxR, int bias, int N, int D, int K, double alpha, double s2Y, double S2B, double d2u, int maxK, int Nsim  )
{
return List::create(B_out, Z_out,theta_out);
}


// initialisation
//RcppGSL::matrix_view Xview = gsl_matrix_view_array(X_dou, D,N);


