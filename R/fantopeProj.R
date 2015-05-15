#' A function that projects a matrix to Fantope
#' @param S An input p by p symmetric matrix
#' @param ndim The dimension of Fantope
#' @return projH The projection matrix of input Smat

fantopeProj <- function(S, ndim){
  # SVD -- need to speed up
  Seig = eigen(S)
  D = Seig$values
  V = Seig$vectors

  # threshold eigen values
  theta = getTheta(D, ndim)$theta
  newD = pmin(pmax(D-theta, 0), 1)

  # reconstruct
  projH = V %*% diag(newD) %*% t(V)
  return(projH)
}
