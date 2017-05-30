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
#include "../core/InferenceFunctions.cpp"

