#' Brief Description
#'
#' Long Description
#'
#'
#' @param x
#' @param y
#'
#' @return double
#' @export
#' @examples
#'

linear_reg <- function(x, y) {
  n <- dim(x)[1]
  p <- dim(x)[2] + 1
  y <- as.matrix(y)
  x1 <- as.matrix(cbind(Intercept = rep(1, n), x))
  xtx_inv <- solve(t(x1) %*% x1)
  coeffs <- xtx_inv %*% t(x1) %*% y
  fv <- x1 %*% coeffs
  res <- y - fv
  s2 <- sum(res^2) / (n - p)
  return(list(coefficients = coeffs, fitted_values = fv, residuals = res,
              sigma_2_hat = s2, xtx_inverse = xtx_inv, n = n, p = p))
}


