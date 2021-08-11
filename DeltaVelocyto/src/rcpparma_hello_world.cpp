// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


//' get the PCA of a matrix
//'
//' @param x matrix
//' @export
// [[Rcpp::export]]
arma::mat getPCA(arma::mat matrix) {

    arma::mat coeff;
    arma::mat score;
    arma::vec latent;
    arma::vec tsquared;

    princomp(coeff, score, latent, tsquared, matrix);

    return score;
}
