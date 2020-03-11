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

lr_coefficient_CI <- function(lr, alpha) {
  qt <- qt(1 - (alpha / 2), (lr$n - lr$p))
  var_coeffs <- lr$sigma_2_hat * diag(lr$xtx_inverse)
  lower <- lr$coefficients - (qt * sqrt(var_coeffs))
  upper <- lr$coefficients + (qt * sqrt(var_coeffs))
  CIs <- cbind(lower, upper)
  colnames(CIs) <- c("Lower", "Upper")
  return(CIs)
}

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
