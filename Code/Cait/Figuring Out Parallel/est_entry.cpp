// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double est_entry(arma::cx_mat V_star_mat, arma::cx_mat V_mat, arma::mat R_mat, int K){

    double x = arma::norm2est(V_star_mat*R_mat*V_mat, 1e-6, 100);

    double output = (1.0/K)*x;

    return output;
}

