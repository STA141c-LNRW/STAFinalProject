#' Brief Description
#'
#' Long Description
#'
#'
#' @param x matrix
#' @param y column vector
#'
#' @return lr, linear model
#' @export
#' @examples
#'

linear_reg <- function(x, y) {
  #observation count
  n <- dim(x)[1]
  #parameter count
  p <- dim(x)[2] + 1

  #reponse variable
  y <- as.matrix(y)

  #explanatory variables
  x1 <- as.matrix(cbind(Intercept = rep(1, n), x))

  # (X'X)^-1
  xtx_inv <- solve(t(x1) %*% x1)

  #Betas
  coeffs <- xtx_inv %*% t(x1) %*% y

  #fitted values
  fv <- x1 %*% coeffs

  #residuals
  res <- y - fv

  #covariance matrix?
  s2 <- sum(res^2) / (n - p)

  #returns an fitted model object, with features: coeffs, fitted values, residuals,sigma hat2, n, p
  return(list(coefficients = coeffs, fitted_values = fv, residuals = res,
              sigma_2_hat = s2, xtx_inverse = xtx_inv, n = n, p = p))
}


