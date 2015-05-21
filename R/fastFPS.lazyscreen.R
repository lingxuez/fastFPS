#' The main function that performs fast FPS on a sequence of lambdas
#' @param S Input symmetric matrix
#' @param ndim The dimension of Fantope
#' @param lambda The smoothing parameter, from large to small
#' @param maxiter The maximum interations
#' @param eps Accuracy for the stopping criterion
#' @param verbose how many outputs wanted
#' @return A list that contains the following objects:

fastFPS.lazyscreen <- function(S, ndim, lambda, maxiter=100, eps=1e-3, verbose=0){
  p <- nrow(S)
  nsol <- length(lambda)

  ## screening
  act_indices <- findActive(S, ndim, min(lambda))
  Sworking <- S[act_indices, act_indices, drop=FALSE]
  p_act = nrow(Sworking)

  if (verbose > 1){
    print(paste(length(act_indices), "active variables"))
  }

  ## initialization
  solutions <- list(ndim=ndim, lambda=lambda,
                    projection=vector("list",nsol),
                    leverage=matrix(0, nrow=p, ncol=nsol),
                    L1=rep(0, nsol),
                    var.explained=rep(0, nsol))
  y0 <- matrix(0, p_act, p_act) # aux variable
  w0 <- matrix(0, p_act, p_act) # dual variable
  tau <- max(abs(Sworking)) # admm penalty


  ## solution path
  for (i in 1:nsol){
    if (verbose > 0){
      cat(".")
    }

    # admm
    sol_admm <- fastADMM(Sworking, ndim, lambda[i], y0, w0, tau,
                   maxiter=maxiter, eps=eps)
    y0 <- sol_admm$y1
    w0 <- sol_admm$w1
    tau <- sol_admm$tau

    # store solutions
    projH = matrix(0, p, p)
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
