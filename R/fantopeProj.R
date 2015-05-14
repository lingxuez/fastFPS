#' A function that projects a matrix to Fantope
#' @param Smat An input p by p symmetric matrix
#' @param ndim The dimension of Fantope
#' @return projH The projection matrix of input Smat

fantopeProj <- function(Smat, ndim){
  # SVD -- need to speed up
  Seig = eigen(Smat)
  D = Seig$values
  V = Seig$vectors

  # threshold eigen values
  theta = getTheta(D, ndim)
  newD = pmin(pmax(D-theta, 0), 1)

  # reconstruct
  projH = V %*% diag(newD) %*% t(V)
  return(projH)
}
