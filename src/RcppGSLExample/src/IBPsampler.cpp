// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-
//
// Copyright (C)  2010 - 2015 Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppGSL.
//
// RcppGSL is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppGSL is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppGSL.  If not, see <http://www.gnu.org/licenses/>.


#include <RcppGSL.h>
#include <Rcpp.h>
#include "GeneralFunctions.h"
#include "InferenceFunctions.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"


// [[Rcpp::export]]
Rcpp::List IBPsampler(const RcppGSL::Matrix input_X, Rcpp::CharacterVector  input_C, RcppGSL::Matrix & input_Z, SEXP input_bias, Rcpp::NumericVector  input_F, SEXP input_s2u, SEXP input_s2B, SEXP input_alpha, SEXP input_Nsim, SEXP input_maxK, SEXP input_missing) {
  
    //..................CHECKING INPUTS AND OUTPUTS.............//
    int D = input_X.nrow(), N = input_X.ncol(); //Input_X should be a matrix of size DxN
    int K = input_Z.nrow(); //Input_Z should be a matrix of size KinixN
    int Cdim= input_C.size();
    int bias = Rcpp::as<int>(input_bias);
    double s2u = Rcpp::as<double>(input_s2u);
    double s2B = Rcpp::as<double>(input_s2B);
    double alpha = Rcpp::as<double>(input_alpha);
    int Nsim = Rcpp::as<int>(input_Nsim);
    int maxK = Rcpp::as<int>(input_maxK);
    double missing = Rcpp::as<double>(input_missing);
    
    if (Cdim!=D){
      throw std::range_error("Invalid number of dimensions for vector C\n");
    }
    if (input_Z.ncol()!=N){
      throw std::range_error("Invalid number of dimensions for initial matrix Z\n");
    }
    
    char C[D];// = Rcpp::as<string>(input_C);
    double  f[D];
    printf("Cvec=  ");
    for (int dd = 0; dd < D; dd++) {
      C[dd]=*(char*)input_C[dd];
      f[dd]=(input_F(dd));
      C[dd] = std::tolower(C[dd]);
      //printf("%c ", *(char*)input_C[dd]);
      printf("%c ", C[dd]);
    }
//    printf("\n N=%d, D=%d, Kini=%d, Cdim=%d \n", N, D, K, Cdim);
    printf("\n ");
    //...............BODY C CODE.......................//
    //Initialization
    gsl_matrix *X = gsl_matrix_alloc(D,N);
    gsl_matrix_memcpy(X , input_X);
    gsl_matrix *Z = gsl_matrix_calloc(maxK,N);
    for (int i=0;i<N;i++){
      for (int k=0; k<K; k++){
        gsl_matrix_set (Z, k, i,gsl_matrix_get (input_Z, k, i));
        //printf("znk=%f",gsl_matrix_get (input_Z, k, i));
      }
    }
    gsl_matrix **B=(gsl_matrix **) calloc(D,sizeof(gsl_matrix*));
    gsl_vector **theta=(gsl_vector **) calloc(D,sizeof(gsl_vector*));
    double w[D],mu[D],s2Y[D];
    int R[D];
    printf("In C++: Transforming input data... ");
    int maxR=initialize_func (N,  D,  maxK, missing,  X, C, B, theta, R, f, mu,  w, s2Y);
    printf("maxR=%d\n", maxR);
    //...............Inference Function.......................//
    printf("In C++: Running Inference Routine... ");
    //int Kest = K;
    int Kest = IBPsampler_func (missing, X, C, Z, B, theta, R, f, mu, w, maxR, bias, N, D, K, alpha, s2B, s2Y, s2u, maxK, Nsim);
    printf("done\n");
    
    //...............SET OUTPUT.......................//
    printf("\n N=%d, D=%d, Kest=%d \n", N, D, Kest);
    
    Rcpp::NumericMatrix out_Z(Kest,N);
    for (int j = 0; j < N; j++) {
        for (int k=0; k<Kest; k++){
          out_Z[Kest*j+k]=gsl_matrix_get (input_Z, k, j);
        }
    }
    
    Rcpp::ListMatrix out_B(D,1);
    int idx_tmp;
    for (int d=0; d<D; d++){
      if (C[d] =='c') {
        idx_tmp = R[d];
       } else {
         idx_tmp = 1;
       }
       Rcpp::NumericMatrix auxB(Kest,idx_tmp); 
       for (int j = 0; j < idx_tmp; j++) {
         for (int k=0; k<Kest; k++){
           auxB[Kest*j+k]=gsl_matrix_get (B[d], k, j);
         }
       }
       out_B[d]=auxB;
    }
    
    Rcpp::NumericMatrix out_theta(D,maxR);      
    for (int d = 0; d < D; d++) {
      for (int i=0;i<maxR;i++){
        if (C[d]=='o' & i<(R[d]-1)){
          out_theta[D*i+d]=gsl_vector_get (theta[d], i);
        }
      }
    }
    
    Rcpp::NumericVector output_MU(D);
    Rcpp::NumericVector output_W(D);
    Rcpp::NumericVector output_s2Y(D);
    
    for (int d=0; d<D; d++){
      output_MU[d]=mu[d];
      output_W[d]=w[d];
      output_s2Y[d]=s2Y[d];
    }
    
    //..... Free memory.....//
    for (int d=0; d<D; d++){
      gsl_matrix_free(B[d]);
      if (C[d] == 'o') { // TODO: Verify why this line gives segmentation fault
        gsl_vector_free(theta[d]);
      }
    }
    gsl_matrix_free(Z);
    free(B);
    free(theta);
    // We return a list with all the outputs
    return Rcpp::List::create(Rcpp::Named("Z") = out_Z, 
                              Rcpp::Named("B")       = out_B,
                              Rcpp::Named("theta")  = out_theta,
                              Rcpp::Named("mu") = output_MU,   
                              Rcpp::Named("w") = output_W, 
                              Rcpp::Named("s2y") = output_s2Y);
}



