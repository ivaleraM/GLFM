#ifndef _INCL_MVNORM_
#define _INCL_MVNORM_
#include <RcppEigen.h>
using namespace Rcpp;

// generate samples from a multivariate normal distribution
NumericMatrix rmvnorm(int nsamples, const NumericVector & mean, const NumericMatrix & covariance);

NumericMatrix centered_rmvnorm(int nsamples, const Eigen::MatrixXd & cholesky_covariance);

// evaluate probability density function of a multivariate normal distribution
NumericVector dmvnorm(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance);

#endif

