#' A function that projects a matrix to Fantope -- speed up version
#' @param S An input p by p symmetric matrix
#' @param ndim The dimension of Fantope
#' @return projH The projection matrix of input Smat

fastFantopeProj <- function(S, ndim){
  p = nrow(S)
  # Zero elements in S?

  # SVD: only values
  Seig = eigen(S, symmetric=TRUE, only.values=TRUE)

  # thres eigen values
  D = Seig$values
  theta = getTheta(D, ndim)
  nActive = theta$last_act_i
  newD = pmin(pmax(D[1:nActive] - theta$theta, 0), 1)

  # eigen vectors
  if (nActive < p){
    V = eigs_sym(S, k=nActive, symmetric=TRUE)$vectors
  } else{
    V = eigen(S, symmetric=TRUE)$vectors
  }

  # reconstruct
  # system.time({
  # projH = V %*% diag(newD) %*% t(V)
  # })
  # system.time({
    projH = reconSVDcpp(diag(newD), V)
  # })
  # projH = reconSVD(diag(newD), V)

  # return(projH)
}
