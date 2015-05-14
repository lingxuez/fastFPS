#' Performs soft-threshold
#' @param x A vector or a matrix
#' @param lambda The threshold
#' @return xThres, the thresholded x

softThreshold <- function(x, lambda){
  xThres = sign(x) * pmax(abs(x)-lambda, 0)
  return(xThres)
}
