
#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double tc(double quant, int tf){
  boost::math::students_t dist(tf);
  return quantile(dist, quant);
}


// [[Rcpp::export]]
NumericVector lr_prediction_interval_C(List lr, arma::colvec x, double alpha){
  x.insert_rows(0, 1);
  x[0] = 1;
  double sigma = lr["sigma_2_hat"];
  NumericMatrix tempinv = lr["xtx_inverse"];
  arma::mat xtx_inverse = arma::mat(tempinv.begin(), tempinv.nrow(),
                                    tempinv.ncol(),
                                    false);
  double var = (1 + trans(x) * xtx_inverse * x)[0] * sigma;
  int n = lr["n"];
  int p = lr["p"];
  double qt = tc(1 - alpha/2, n - p);
  NumericMatrix temp = lr["coefficients"];
  double mean_res = 0;
  for (unsigned int i = 0; i < x.n_rows; i++){
    mean_res += temp[i]*x[i];
  }
  double lower = mean_res - (qt * sqrt(var));
  double upper = mean_res + (qt * sqrt(var));
  NumericVector final = NumericVector::create(lower,mean_res,upper);
  final.attr("dim") = Dimension(1,3);
  return final;
}

// [[Rcpp::export]]
NumericVector lr_prediction_interval_C2(List lr, NumericVector x, double alpha){
  x.insert(0, 1);
  double sigma = lr["sigma_2_hat"];
  NumericMatrix tempinv = lr["xtx_inverse"];
  NumericVector y(tempinv.ncol());
  for (int i = 0; i < tempinv.ncol(); i++){
    y[i] = std::inner_product(x.begin(), x.end(),
                       tempinv.begin()+i*tempinv.ncol(), 0.);
  }
  double t = std::inner_product(y.begin(), y.end(),
                                   x.begin(), 0.);
  double var = (1 + t) * sigma;
  int n = lr["n"];
  int p = lr["p"];
  double qt = tc(1 - alpha/2, n - p);
  NumericMatrix temp = lr["coefficients"];
  double mean_res = 0;
  for (unsigned int i = 0; i < x.size(); i++){
    mean_res += temp[i]*x[i];
  }
  double lower = mean_res - (qt * sqrt(var));
  double upper = mean_res + (qt * sqrt(var));
  NumericVector final = NumericVector::create(lower,mean_res,upper);
  final.attr("dim") = Dimension(1,3);
  return final;
}



