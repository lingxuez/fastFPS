## Check that fast FPS outputs same results as Vu's fps

test_that("fastFPS() works for sample covariance matrix of X~N(0, I_p),
          p=30, n=100, ndim=3",{
  p = 100
  set.seed(100)
  X = matrix(rnorm(p*10), nrow=10, ncol=p)
  S = cov(X)
  ndim = 3
  lambda = quantile(abs(S), probs=c(0.95, 0.9, 0.85, 0.8))
  vu <- fps::fps(S, ndim=ndim, lambda=lambda)
  my <- fastFPS(S, ndim=ndim, lambda=lambda)
  my.lazy <- fastFPS.lazyscreen(S, ndim=ndim, lambda=lambda)
  my.cpp <- fastFPS_cpp(S, ndim=ndim, lambda=lambda)

  all.equal(vu$projection, my$projection, tolerance=0.01)
  all.equal(vu$projection, my.lazy$projection, tolerance=0.01)
  all.equal(vu$projection, my.cpp$projection, tolerance=0.01)
})
