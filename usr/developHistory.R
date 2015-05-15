library("devtools")
library("roxygen2")

# dependencies
devtools::use_package("fps")
devtools::use_package("ggplot2")

# unit test
devtools::use_testthat()
devtools::test()

# develop
devtools::document()
load_all()
test()

# compare with fps package
library(fps)
rm(list=ls())
set.seed(100)
X = matrix(rnorm(900), nrow=10, ncol=90)
S = cov(X)
ndim = 3
lambdas = c(1:5)*0.1
system.time({
  myResult = fastFPS(S, ndim=ndim, lambdas=lambdas, verbose=3)
})

system.time({
  trueResult = fps(S, ndim=ndim, lambda=lambdas, verbose=3)
})

plot(myResult$projection[[3]], trueResult$projection[[3]])
plot(myResult$leverage, trueResult$leverage)

abline(0,1 )
vu = fps(S, ndim=ndim, lambda=lambda, verbose=3)


Rprof("./usr/fastFPSprof.out")
myResult = fastFPS(S, ndim=ndim, lambdas=lambdas, verbose=3)
Rprof(NULL)
summaryRprof("./usr/fastFPSprof.out")
