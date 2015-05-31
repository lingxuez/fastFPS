#' A function that mimics Vu's simplex_sum function
#'
#' Calculated the thresholded sum of a vector
#' @param v A p by 1 vector of sorted eigen values, v_1 >= ... >= v_p
#' @param theta The threshold
#' @return v_sum, the thresholded sum of vector v^+,
#' where v^+ = min(max(x-theta, 0), 1)

simplex_sum <- function(v, theta){
  v = v-theta
  v[v>1] <- 1
  v[v<0] <- 0
  return (sum(v))
}

#' A function that mimics Vu's simplex function
#'
#' Should be working the same way as getTheta
#' @param v A p by 1 vector of sorted eigen values, v_1 >= ... >= v_p
#' @param ndim The dimension of simplex
#' @return theta, the threshold
#' where v^+ = min(max(x-theta, 0), 1)

simplex <- function(v, ndim){
  knots = sort(unique(c(v-1, v)), decreasing=TRUE)

  for (dd in 1:length(knots) ){
    f = simplex_sum(v, theta=knots[dd])
    if (f >= ndim){
      break
    }
  }

  # interpolate; true theta is between [ knots[dd], knots[dd-1] )
  lower = knots[dd];
  upper = knots[dd-1];
  f.lower = simplex_sum(v, lower);
  f.upper = simplex_sum(v, upper);
  theta = lower + (upper-lower) * (ndim-f.lower) / (f.upper-f.lower);

  # rank
  rank = length(v) - sum(v <= theta)

  return (list(theta=theta, rank=rank))

}
