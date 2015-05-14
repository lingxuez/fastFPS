#' Performs ADMM to solve FPS
#' @param S Input symmetric matrix
#' @param ndim The dimension of Fantope
#' @param lambda The smoothing parameter
#' @param maxiter The maximum interations
#' @param eps Accuracy for the stopping criterion
#' @return H The global solution matrix H for the Fantope problem

admm <- function(S, ndim, lambda, maxiter=100, eps=10^(-6)){
  p = nrow(S)
  tau = 0.1* max(abs(S)) # admm penalty
  tauStep = 2 # step size to adjust tau

  y0 = matrix(0, nrow=p, ncol=p)
  w0 = matrix(0, nrow=p, ncol=p)

  niter = 0
  maxnorm = eps+1
  while (niter < maxiter && maxnorm > eps){
    # update
    H = fantopeProj(y0-w0+S/tau, ndim)
    y1 = softThreshold(H+w0, lambda/tau)
    w1 = w0 + H - y1

    # stopping criterion
    normr1 = (Matrix::norm(H-y1, "F"))^2
    norms1 = (Matrix::norm(tau*(y0-y1), "F"))^2
    maxnorm = max(normr1, norms1)
    niter = niter+1

    # update
    y0 = y1
    if (normr1>100*norms1){
      tau = tau*tauStep
      w0 = w1/tauStep
    } else if (norms1 > 100*normr1){
      tau = tau/tauStep
      w0 = w1*tauStep
    } else{
      w0 = w1
    }
  }

  return (y1)
}


