#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/

//[[Rcpp::export]]
NumericVector sammy(int x, NumericVector y){
  NumericVector sample;
  for(int i=0; i<x; i++){
    int index=rand() % y.size();
    sample.insert(i,y[index]);
  }
  return(sample);
}

/***R
sammy(3,10:20)
*/








