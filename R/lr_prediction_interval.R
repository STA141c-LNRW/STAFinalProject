#' Brief Description
#'
#' Long Description
#'
#'
#' @param lr linear model
#' @param x
#' @param alpha significance level
#'
#' @return double
#' @export
#' @examples
#'
#'
lr_prediction_interval <- function(lr, x, alpha) {
  x <- c(1, x)
  var <- lr$sigma_2_hat * (1 + (t(x) %*% lr$xtx_inverse %*% x))
  qt <- qt(1 - (alpha / 2), (lr$n - lr$p))
  mean_res <- x %*% lr$coefficients
  lower <- mean_res - (qt * sqrt(var))
  upper <- mean_res + (qt * sqrt(var))
  PI <- cbind(lower, mean_res, upper)
  colnames(PI) <- c("Lower", "Fit", "Upper")
  return(PI)
}
