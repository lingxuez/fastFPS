#' Computes a set of indices that contains the active variables of FPS
#'
#' If maxoffdiag(i) <= lambda _and_ |{j: diag(j) >= diag(i)}| >= ndim,
#' then variable i can be safely eliminated.
#' @param S The input symmetric matrix S
#' @param ndim The dimension of Fantole
#' @param lambda The smoothing parameter
#' @return act_indices, the indices of active set

findActivePrep <- function(S, ndim){
  p = nrow(S)

  Sdiag <- diag(S)
  S <- abs(S)
  diag(S) <- NA
  Smaxoff = apply(S, 1, max, na.rm=TRUE)

  # sort (diag, maxoffdiag) in lexicographical and descending order
  Spair = cbind(index=c(1:p), Sdiag, Smaxoff)
  lex.order <- order(Sdiag, Smaxoff, decreasing=TRUE)
  Spair <- Spair[lex.order, ]

  return(Spair)
}

