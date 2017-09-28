#include <Rcpp.h>
using namespace Rcpp;

extern "C"
void mcCuda(double *vec, double *samps, double *mat, int N, int S);

//' Monte Carlo estimated covariance, CUDA accelerated
//'
//' @param x a vector of distances
//' @param s a vector of random samples from the estimated spectral density
//' @return a vector or some shit I dunno
// [[Rcpp::export]]
NumericVector mc(NumericVector x, NumericVector s) {
  int n = x.size();
  int m = s.size();
  NumericVector result = NumericVector(n);
  double *xx = x.begin();
  double *ss = s.begin();
  double *rr = result.begin();
  mcCuda(xx, ss, rr, n, m);
  return result;
}
