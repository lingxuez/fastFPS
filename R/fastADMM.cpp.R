#' Performs ADMM to solve FPS, using fast fantope projection
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

fastADMM.cpp <- function(S, ndim, lambda, y, w, tau,
                     tauStep=2, maxiter=100, eps=10^(-3)){
  eps = sqrt(ndim)*eps # stopping criterion

  for (niter in 1:maxiter){
    y_old = y

    # update
    H = fastFantopeProj(y-w+S/tau, ndim)
    y = softThresholdCpp(H+w, lambda/tau)
    w = w + H - y

    # stopping criterion
    normr1 = Matrix::norm(H-y, "F")
    norms1 = tau * Matrix::norm(y-y_old, "F")
    if (normr1 < eps & norms1 < eps){
      break
    }

    # update
    if (normr1 > 10*norms1){
      tau = tau*tauStep
      w = w/tauStep
    } else if (norms1 > 10*normr1){
      tau = tau/tauStep
      w = w*tauStep
    }

  }

  return(list(y=y, w=w, tau=tau, niter=niter))
}


