// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double KS(arma::colvec x, arma::colvec y) {
  int n = 100;
  arma::colvec cc = join_cols(x, y);
  double minx = cc.min();
  double maxx = cc.max();
  arma::colvec r = arma::linspace(minx, maxx, n);
  int N = r.n_rows;

  arma::colvec fx = arma::linspace(minx, maxx, n);
  arma::colvec v = arma::sort(x);
  arma::colvec s = v;
    for (int i=0; i<n;i++) {
      double X = r(i);
      s.fill(0.0); s.elem( find(v < X) ).ones();
      fx(i) = mean(s);
    }

  arma::colvec fx1 = arma::linspace(minx, maxx, n);
  arma::colvec v1 = arma::sort(y);
  arma::colvec s1 = v1;
    for (int i=0; i<n;i++) {
      double X = r(i);
      s1.fill(0.0); s1.elem( find(v1 < X) ).ones();
      fx1(i) = mean(s1);
    }

    arma::colvec diff = arma::abs(fx-fx1);
  return arma::max(diff);
}
// [[Rcpp::export]]
Rcpp::NumericVector K_S(arma::mat mt,arma::mat mt2) {
  int n = mt.n_cols;
  Rcpp::NumericVector results(n);
  for (int i=0; i<n;i++) {
    arma::colvec x=mt.col(i);
    arma::colvec y=mt2.col(i);
    results[i] = KS(x, y);
  }
  return results;
}


