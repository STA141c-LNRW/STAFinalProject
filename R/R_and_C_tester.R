library(Rcpp)
library(bench)
sourceCpp("c_script.cpp") #may use this file in future when implementing c++ functions



testdata1 <- rnorm(1e5,50,10)

bench::mark(sample(testdata1, 100))
bench::mark(sammy(testdata1, 100))
#Rcpp version uses less memory and is faster
#For larger sample sizes, the standard R version starts to become faster and more memory efficient

