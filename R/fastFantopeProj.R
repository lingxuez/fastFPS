#' A function that projects a matrix to Fantope -- speed up version
#' @param S An input p by p symmetric matrix
#' @param ndim The dimension of Fantope
#' @return projH The projection matrix of input Smat

fastFantopeProj <- function(S, ndim){
  p = nrow(S)

  # SVD: only values
  D = eigen(S, symmetric=TRUE, only.values=TRUE)$values

  # thres eigen values
  theta = getTheta(D, ndim)
  nActive = theta$last_act_i

  # eigen vectors
  if (nActive < p & p >= 3){
    V = rARPACK::eigs_sym(S, k=nActive, symmetric=TRUE)$vectors
    newD = pmin(pmax(D[1:nActive] - theta$theta, 0), 1)
  } else{
    V = eigen(S, symmetric=TRUE)$vectors
    newD = D
  }

  # reconstruct
  if (p==1){
    projH = V * newD * V
  } else {
    projH = V %*% diag(newD) %*% t(V)
  }

  return(projH)
}
