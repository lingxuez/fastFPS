#' The main function that performs fast FPS on a sequence of lambdas
#' @param S Input symmetric matrix
#' @param ndim The dimension of Fantope
#' @param lambda The smoothing parameter
#' @param maxiter The maximum interations
#' @param eps Accuracy for the stopping criterion
#' @return A list that contains the following objects:

fastFPS <- function(S, ndim, lambdas, maxiter=100, eps=1e-3, verbose=0){
  p <- nrow(S)
  nsol <- length(lambdas)

  ## screening
  act_indices <- findActive(S, ndim, min(lambdas))
  Sworking <- S[act_indices, act_indices]
  p_act <- nrow(Sworking)
  if (verbose > 0){
    print(paste("After screening,", p_act, "variables are remained active"))
  }

  ## initialization
  solutions <- list(ndim=ndim, lambda=lambdas,
                   projection=vector("list",nsol),
                   leverage=matrix(0, nrow=p, ncol=nsol),
                   L1=rep(0, nsol),
                   var.explained=rep(0, nsol))
  y0 <- matrix(0, p_act, p_act) # aux variable
  w0 <- matrix(0, p_act, p_act) # dual variable
  tau <- max(abs(Sworking)) # admm penalty
  projH <- matrix(0, p, p)

  ## solution path
  for (i in 1:nsol){
    if (verbose > 0){
      cat (".")
    }

    # admm
    sol_admm <- admm(Sworking, ndim, lambdas[i], y0, w0, tau)
    y0 <- sol_admm$y1
    w0 <- sol_admm$w1
    tau <- sol_admm$tau

    # store solutions
    projH = matrix(0, p, p)
    attributes(projH)$dimnames = attributes(S)$dimnames
    projH[act_indices, act_indices] = y0
    solutions$projection[[i]] = projH

    solutions$leverage[,i] = diag(projH)
    solutions$L1[i] = sum(rowSums(abs(sol_admm$y1))) # |H|_1
    solutions$var.explained[i] = sum(rowSums(Sworking * y0)) # <S, H>

    if (verbose > 1){
      cat(sol_admm$niter)
    }
    if (verbose > 2){
      cat(paste0("(", round(sol_admm$tau, 5), ")"))
    }
  }

  return(solutions)
}
