#' Performs ADMM to solve FPS
#'
#' Solve: min_{H, y} infty*I_{Fd}(H) - <S,H> + lambda*|y|_1 + tau/2 |H-y|_F^2
#' subject to H-y=0
#' @param S Input symmetric matrix
#' @param ndim The dimension of Fantope
#' @param lambda The smoothing parameter
#' @param maxiter The maximum interations
#' @param eps Accuracy for the stopping criterion
#' @param y0 Initialized value for aux variable
#' @param w0 Initialized value for dual variable
#' @param tau ADMM penalty
#' @param tauStep Step size to adjust tau
#' @return H The global solution matrix H for the Fantope problem

admm <- function(S, ndim, lambda, y0, w0, tau, tauStep=2, maxiter=100, eps=10^(-3)){
  eps = sqrt(ndim)*eps # stopping criterion

  niter = 0
  maxnorm = eps+1
  while (niter < maxiter && maxnorm > eps){
    # update
    H = fantopeProj(y0-w0+S/tau, ndim)
    y1 = softThreshold(H+w0, lambda/tau)
    w1 = w0 + H - y1

    # stopping criterion
    normr1 = Matrix::norm(H-y1, "F")
    norms1 = tau * Matrix::norm((y0-y1), "F")
    maxnorm = max(normr1, norms1)
    niter = niter+1

    # update
    y0 = y1
    if (normr1 > 10*norms1){
      tau = tau*tauStep
      w0 = w1/tauStep
    } else if (norms1 > 10*normr1){
      tau = tau/tauStep
      w0 = w1*tauStep
    } else{
      w0 = w1
    }
  }

  return(list(y1=y1, w1=w1, tau=tau, niter=niter))
}


