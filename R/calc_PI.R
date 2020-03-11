#' Brief Description
#'
#' Long Description
#'
#'
#' @param df dataframe
#' @param x_h double float
#'
#' @return interval
#' @export
#' @examples
#'

#include in package
library(tibble)
library(car)
library(purrr)
library(apricom)


#Note: it functions more like a confidence interval for predictions
#USE: predict(fit,data.frame(hp = 10, disp = 10), interval = "confidence")
#as a demonstration
calc_PI = function(df,x_1,x_2){

  #calculates individual forecast
  calc_forecast = function(df, x_1,x_2){
    fit = lm(mpg ~ hp + disp, data = df)
    beta_0 = fit$coefficients[[1]]
    beta_1 = fit$coefficients[[2]]
    beta_2 = fit$coefficients[[3]]
    (y_hat = beta_0 + beta_1*x_1 + beta_2*x_2)
  }
  #calculates estimates of forecasts using bootstrapping
  boot_forecast = function(i, df, x_1,x_2){
    index = sample(x = seq_len(100), size = 100, replace = TRUE)
    sub_data = df[index,]
    calc_forecast(sub_data, x_1, x_2)
  }

  r = 100
  #apply x_h
  #prediction interval
  pi_list = map_dbl(seq_len(r), boot_forecast, df, x_1, x_2)
  #note data assumptions

  quantile(pi_list,c(0.025,0.975))
}
