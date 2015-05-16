#' Computes a set of indices that contains the active variables of FPS
#'
#' If maxoffdiag(i) <= lambda _and_ |{j: diag(j) >= diag(i)}| >= ndim,
#' then variable i can be safely eliminated.
#' @param S The input symmetric matrix S
#' @param ndim The dimension of Fantole
#' @param lambda The smoothing parameter
#' @return act_indices, the indices of active set

findActive <- function(S, ndim, lambda){
  p = nrow(S)
  Sdiag <- diag(S)
  Soff <- S
  diag(Soff) <- NA
  Smaxoff = apply(Soff, 1, max, na.rm=TRUE)

  # sort (diag, maxoffdiag) in lexicographical and descending order
  Spair = cbind(index=c(1:p), Sdiag, Smaxoff)
  lex.order <- order(Sdiag, Smaxoff, decreasing=TRUE)
  Spair <- Spair[lex.order, ]
  # active set contains those with: j<ndim or maxoff > lambda
  act_rows <- which(c(1:p) < ndim | Spair[, "Smaxoff"] > lambda)
  act_indices = Spair[act_rows, "index"]

  return(act_indices)
}

