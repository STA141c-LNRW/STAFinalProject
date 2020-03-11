#' Brief Description
#'
#' Long Description
#'
#'
#' @param df dataframe
#'
#' @return double
#' @export
#' @examples
#'

calc_slope <- function(df) {
  #explanatory and response variables
  x = df$x
  y = df$y

  #means
  x_bar = mean(x)
  y_bar = mean(y)

  #sum of squares
  sxy = sum((x-x_bar)*(y-y_bar))
  sxx = sum((x-x_bar)^2)

  #beta calculation
  (beta = sxy/sxx)
}
