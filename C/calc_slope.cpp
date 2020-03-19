

#include <Rcpp.h>

using namespace Rcpp;



// [[Rcpp::export]]
double calc_slopeC(DataFrame df){
  NumericVector x = df["x"];
  int length = x.size();
  NumericVector y = df["y"];
  double meanX = 0;
  double meanY = 0;
  for (int i = 0; i < length; i++){
    meanX += x[i];
    meanY += y[i];
  }
  meanX /= length;
  meanY /= length;
  double sumSXX = 0;
  double sumSXY = 0;
  for (int i = 0; i < length; i++){
    double temp = x[i] - meanX;
    sumSXX += temp*temp;
    sumSXY += temp*(y[i] - meanY);
  }
  return sumSXY/sumSXX;
}