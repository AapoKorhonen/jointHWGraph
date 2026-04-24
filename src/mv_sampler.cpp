#include <RcppArmadillo.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

// Function to sample from a multivariate normal distribution
// [[Rcpp::export]]
arma::mat mvrnorm_cpp(int n, const arma::vec& mu, const arma::mat& Sigma) {
  int k = mu.size();
  arma::mat Y = arma::randn(n, k);
  arma::mat L = arma::chol(Sigma, "lower");
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}
