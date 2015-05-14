#' A function that calculates the threshold theta for Fantope projection
#' @param v A p by 1 vector of sorted eigen values, v_1 >= ... >= v_p
#' @param ndim The effective dimension desired
#' @return theta, a threshold such that /sum_i v_i^+ = ndim/,
#' where x^+ = min(max(x-theta, 0), 1)

getTheta <- function(v, ndim){

  # if /v_{ndim} - v_{ndim+1} >= 1/, then simply let theta=v_{ndim}-1,
  # will have v_1^+ = ... = v_{ndim}^+ = 1, v_{ndim+1}^+ = ... = v_p^+ = 0
  if (v[ndim]-v[ndim+1] >= 1){
    return (v[ndim]-1)
  }

  p = length(v)
  v1 = c(v, v[p]-ndim/p) # v_p - ndim/p is a lowerbound for theta

  # Search whether theta=v_i for some i
  ddnew = 0
  dnew = max(ndim-2, 0) # start from v_{ndim-1} or v_1 -- why not ndim?
  fnew = 0
  while (fnew < ndim){
    # store the value from last v_{i-1}
    f = fnew
    d = dnew
    dd = ddnew

    # move to the next v_i
    dnew = d + 1
    theta = v1[dnew]
    ddnew = which(v1-theta < 1)[1]
    # v_i^+ = 0 for i<= dnew,
    # v_i^+ = v_i-theta for dnew < i <= ddnew
    # v_i^+ = 1 for i > ddnew
    fnew = sum(v1[dnew:ddnew] - theta) + (ddnew-1)
  }

  if (fnew == ndim){
    return (theta)
  }

  # else, theta must be between (v1[d+1], v1[d]), where dnew=d+1
  theta = v1[d]
  # step size: each time release v_dd to be 1, while remaining theta>=v1[d+1]
  m0 = min(theta-v1[d+1], 1-(v1[dd]-theta))
  while ((f + (d-dd+1)*m0) < ndim){
    f = f + (d-dd+1)*m0
    theta = theta - m0
    dd = dd + 1
    m0 = min(theta-v1[d+1], 1-(v1[dd]-theta))
  }

  # real theta must be between (theta-m0, theta)
  theta = theta - (ndim-f)/(d-dd+1)
  return (theta)
}
