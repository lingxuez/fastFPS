#' Computes a set of indices that contains the active variables of FPS
#'
#' If maxoffdiag(i) <= lambda _and_ |{j: diag(j) >= diag(i)}| >= ndim,
#' then variable i can be safely eliminated.
#' @param S The input symmetric matrix S
#' @param ndim The dimension of Fantole
#' @param lambda The smoothing parameter
#' @return act_indices, the indices of active set

findActiveSeq <- function(S, ndim, lambda){
  p = nrow(S)
  nsol = length(lambda)

  Sdiag <- diag(S)
  S.abs.off <- abs(S)
  diag(S.abs.off) <- NA
  Smaxoff = apply(S.abs.off, 1, max, na.rm=TRUE)

  # sort (diag, maxoffdiag) in lexicographical and descending order
  Spair = cbind(index=c(1:p), Sdiag, Smaxoff)
  lex.order <- order(Sdiag, Smaxoff, decreasing=TRUE)
  Spair <- Spair[lex.order, ]

  # active set contains those with: j<=ndim or maxoff > lambda
  act_indicesSeq <- vector("list", nsol)
  for (i in 1:nsol){
    act_rows <- which(c(1:p) <= ndim | Spair[, "Smaxoff"] > lambda[i])
    act_indicesSeq[[i]] = Spair[act_rows, "index"]
  }

  return(act_indicesSeq)
}

