library("devtools")
library("roxygen2")

# dependencies
devtools::use_package("fps")
devtools::use_package("ggplot2")
devtools::use_package("Matrix")
devtools::use_package("rARPACK")
devtools::use_package("inline")
devtools::use_package("RcppArmadillo")

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
p = 300
set.seed(100)
X = matrix(rnorm(p*10), nrow=10, ncol=p)
S = cov(X)
ndim = 3
lambdas = sort(quantile(S[lower.tri(S)], probs=c(0.8, 0.85, 0.9, 0.95, 0.99)), decreasing=TRUE)

plot(myResult$projection[[1]], trueResult$projection[[1]])
plot(myResult$leverage, trueResult$leverage)

abline(0,1 )
vu = fps(S, ndim=ndim, lambda=lambda, verbose=3)

Rprof("./usr/myFPSprof.out")
system.time({
  myResult = myFPS(S, ndim=ndim, lambdas=lambdas, verbose=3)
})
Rprof(NULL)
summaryRprof("./usr/myFPSprof.out")

Rprof("./usr/fastFPSprof.out")
system.time({
  myFastResult = fastFPS(S, ndim=ndim, lambdas=lambdas, verbose=3)
})
Rprof(NULL)
summaryRprof("./usr/fastFPSprof.out")

system.time({
  trueResult = fps(S, ndim=ndim, lambda=lambdas, verbose=3)
})

plot(S, myResult$projection[[3]], pch=".")
abline(v=lambdas[3])

H = myResult$projection[[3]]
sum(abs(H[lower.tri(H)]) < 1e-3 ) / (90*89/2)

##############
## Try RcppArmadillo
##############
library(inline)
library(RcppArmadillo)
g <- cxxfunction ( signature (vs = "numeric"),
                   plugin = "RcppArmadillo", body = '
                   arma::vec v = Rcpp::as < arma::vec > (vs);
                   arma::mat op = v * v. t () ;
                   double ip = arma::as_scalar ( v.t () * v);
                   return Rcpp::List::create ( Rcpp::Named ("outer") =op,
                                                  Rcpp::Named ("inner") = ip );
                   ')

g(7:11)




