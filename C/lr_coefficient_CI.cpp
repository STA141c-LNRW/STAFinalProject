
#include <Rcpp.h>
#include <boost/math/distributions/students_t.hpp>



// [[Rcpp::depends(BH)]]

using namespace Rcpp;
// [[Rcpp::export]]
double tc(double quant, int tf){
  boost::math::students_t dist(tf);
  return quantile(dist, quant);
}

// [[Rcpp::export]]
SEXP tr(double quant, int tf){
  Function qt("qt");
  return qt(quant, tf);
}

// [[Rcpp::export]]
NumericVector lr_coefficient_CI_C(List lr, double alpha){
  int n = lr["n"];
  int p = lr["p"];
  double q = tc(1 - alpha/2, n - p);
  NumericMatrix inverse = lr["xtx_inverse"];
  double sigma = lr["sigma_2_hat"];
  NumericVector coefficients = lr["coefficients"];
  NumericMatrix final(coefficients.size(), 2);
  for (int i = 0; i < coefficients.size(); i++){
    double var_coeffs = sigma * inverse(i, i);
    final[i] = coefficients[i] - q*sqrt(var_coeffs);
    final[i+coefficients.size()] = coefficients[i] + q*sqrt(var_coeffs);
  }
  colnames(final) = CharacterVector::create("Lower","Upper");
  return final;
}
