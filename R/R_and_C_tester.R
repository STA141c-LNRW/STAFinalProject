library(Rcpp)
library(bench)
sourceCpp("c_script.cpp")



testdata1 <- rnorm(1e5,50,10)

bench::mark(sample(testdata1, 100))
bench::mark(sammy(testdata1, 100))
#Rcpp version uses less memory and is faster

bench::mark(sample(testdata1, 1000))
bench::mark(sammy(testdata1, 1000))
#For larger sample sizes, the standard R version starts to become faster and more memory efficient

bench::mark(rnorm(1000, 100, 10))
bench::mark(cnorm(1000, 100, 10))
#cnorm appears to be faster

bench::mark(rnorm(4000, 100, 10))
bench::mark(cnorm(4000, 100, 10))
#right around sample sizes of 4000, cnorm starts to take up too much memory and is slower

bench::mark(rgamma(100, 10, 5))
bench::mark(cgamma(100, 10, 5))
# the base r version handles this distribution sampling better than the c++ version


