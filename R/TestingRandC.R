


library(Rcpp)
library(rbenchmark)
set.seed(121)

Random = runif(200000,0,100)
AT = data.frame(x = 1:200000, y = 1:200000 + Random)


x1 = runif(200000,0,1000) + runif(200000, 500, 1000)
x2 = runif(200000,0,2000) + runif(200000, 1000, 2000)
x3 = runif(200000,2000,3000) + runif(200000, 1000, 2000)
x4 = runif(200000,2500,2750) + runif(200000, 1000, 2000)
y = x1 + x2 + x3 + x4 + Random^2 + runif(20000, 500, 777)
xFrame = as.matrix(data.frame(x1,x2,x3,x4))

sourceCpp("R/TestingC.cpp")
source("R/calc_slope.R")



benchmark("C++" = {calc_slopeC(AT)},
          "R" = {calc_slope(AT)},
          replications = 1000,
          columns = c("test","replications","elapsed","relative","user.self","sys.self"))
#p = 4
benchmark(".lm.fit" = {.lm.fit(as.matrix(data.frame(1, xFrame)), y)},
          "lm.fit" = {lm.fit(as.matrix(data.frame(1, xFrame)), y)},
          "C++" = {linear_regC(xFrame, y)},
          "R" = {linear_reg(xFrame, y)},
            replications = 100,
            columns = c("test","replications","elapsed","relative","user.self","sys.self"))

# p =1
benchmark(".lm.fit" = {.lm.fit(as.matrix(data.frame(1, AT[[1]])), as.vector(AT[[2]]))},
          "lm.fit" = {lm.fit(as.matrix(data.frame(1, AT[[1]])), as.vector(AT[[2]]))},
          "C++ with std::inner_product" = {linear_regC(as.matrix(AT[[1]]), as.vector(AT[[2]]))},
          "C++ without std::inner_product" = {linear_regC2(as.matrix(AT[[1]]), as.vector(AT[[2]]))},
          "R" = {linear_reg(as.matrix(AT[[1]]), as.vector(AT[[2]]))},
          replications = 100,
          columns = c("test","replications","elapsed","relative","user.self","sys.self"))

#p = 40
x1 = runif(8000000,0,1000) + runif(8000000, 500, 1000)
FourtyX = matrix(x1, nrow = 200000, ncol = 40)
benchmark(".lm.fit" = {.lm.fit(as.matrix(data.frame(1, FourtyX)), y)},
          "lm.fit" = {lm.fit(as.matrix(data.frame(1, FourtyX)), y)},
          "C++ with std::inner_product" = {linear_regC(FourtyX, y)},
          "C++ without std::inner_product" = {linear_regC2(FourtyX, y)},
          "R" = {linear_reg(FourtyX, y)},
          replications = 1,
          columns = c("test","replications","elapsed","relative","user.self","sys.self"))

benchmark("%*%" = {t(FourtyX) %*% FourtyX},
          "multiply" = {multiply(t(FourtyX), FourtyX)},
          "multiply with std::inner_product" = {multiply2(t(FourtyX), FourtyX)},
          replications = 1,
          columns = c("test","replications","elapsed","relative","user.self","sys.self"))

benchmark("%*%" = {t(xFrame) %*% xFrame},
          "multiply" = {multiply(t(xFrame), xFrame)},
          "multiply with std::inner_product" = {multiply2(t(xFrame), xFrame)},
          replications = 100,
          columns = c("test","replications","elapsed","relative","user.self","sys.self"))

