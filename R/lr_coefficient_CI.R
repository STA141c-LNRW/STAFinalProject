#' Calculates the
#'
#' Long Description
#'
#'
#' @param lr linear model
#' @param alpha significant level
#'
#' @return confidence interval
#' @export
#' @examples
#'
lr_coefficient_CI <- function(lr, alpha) {
  qt <- qt(1 - (alpha / 2), (lr$n - lr$p))
  var_coeffs <- lr$sigma_2_hat * diag(lr$xtx_inverse)
  lower <- lr$coefficients - (qt * sqrt(var_coeffs))
  upper <- lr$coefficients + (qt * sqrt(var_coeffs))
  CIs <- cbind(lower, upper)
  colnames(CIs) <- c("Lower", "Upper")
  return(CIs)
}
