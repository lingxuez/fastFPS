#' The main function that performs FPS on a sequence of lambdas
#'
#' This is an exact replication of Vincent Vu's fps package
#' @param S Input symmetric matrix
#' @param ndim The dimension of Fantope
#' @param lambdas The smoothing parameter, from large to small
#' @param maxiter The maximum interations
#' @param eps Accuracy for the stopping criterion
#' @param verbose How many outputs are wanted
#' @return A list that contains the following objects:

myFPS <- function(S, ndim, lambdas, maxiter=100, eps=1e-3, verbose=0){
  p <- nrow(S)
  nsol <- length(lambdas)

  ## initialization
  solutions <- list(ndim=ndim, lambda=lambdas,
                   projection=vector("list",nsol),
                   leverage=matrix(0, nrow=p, ncol=nsol),
                   L1=rep(0, nsol),
                   var.explained=rep(0, nsol))
  projH <- matrix(0, p, p)
  y0 <- matrix(0, p, p) # aux variable
  w0 <- matrix(0, p, p) # dual variable
  tau <- max(abs(S)) # admm penalty

  ## solution path
  for (i in 1:nsol){
    if (verbose > 0){
      cat (".")
    }

    ## screening
    act_indices <- findActive(S, ndim, lambdas[i])
    Sworking <- S[act_indices, act_indices]
    if (verbose > 0){
      cat("(", length(act_indices), "variables)")
    }

    # admm
    sol_admm <- admm(Sworking, ndim, lambdas[i],
                     y0[act_indices, act_indices], w0[act_indices, act_indices], tau)
    y0[act_indices, act_indices] <- sol_admm$y1
    w0[act_indices, act_indices] <- sol_admm$w1
    tau <- sol_admm$tau

    # store solutions
    projH = matrix(0, p, p)
    attributes(projH)$dimnames = attributes(S)$dimnames
    projH[act_indices, act_indices] = sol_admm$y1
    solutions$projection[[i]] = projH

    solutions$leverage[,i] = diag(projH)
    solutions$L1[i] = sum(rowSums(abs(sol_admm$y1))) # |H|_1
    solutions$var.explained[i] = sum(rowSums(Sworking * sol_admm$y1)) # <S, H>

    if (verbose > 1){
      cat(sol_admm$niter)
    }
    if (verbose > 2){
      cat(paste0("(", round(sol_admm$tau, 5), ")\n"))
    }
  }

  return(solutions)
}
