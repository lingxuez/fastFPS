library("devtools")
library("roxygen2")
library(fps)
library(microbenchmark)

# dependencies
devtools::use_package("fps")
devtools::use_package("ggplot2")
devtools::use_package("Matrix")
devtools::use_package("rARPACK")
devtools::use_package("inline")
devtools::use_package("RcppArmadillo")
devtools::use_rcpp()

# unit test
devtools::use_testthat()
devtools::test()

# develop
devtools::document()
load_all()
test()
check()
# devtools::build(binary=TRUE)
devtools::install("/Users/lingxue/Documents/Thesis/Rtools/fastFPS")

# compare with fps package
rm(list=ls())
p = 100
set.seed(100)
X = matrix(rnorm(p*10), nrow=10, ncol=p)
S = cov(X)
ndim = 3
lambdas = sort(quantile(S[lower.tri(S)], probs=c(0.8, 0.85, 0.9, 0.95, 0.99)), decreasing=TRUE)

checkFpsEqual <- function(values){
  all.equal(values[[1]]$projection, values[[2]]$projection)
}

load_all()
microbenchmark(vu <- fps(S, ndim=ndim, lambda=lambdas, verbose=0),
               my <- fastFPS(S, ndim=ndim, lambda=lambdas, verbose=0),
               check=checkFpsEqual, times=2)

microbenchmark(vu <- fps(S, ndim=ndim, lambda=lambdas[1], verbose=3),
               my <- fastFPS(S, ndim=ndim, lambda=lambdas[1], verbose=3),
               check=checkFpsEqual, times=2)


i.test = 5
all.equal(vu$projection[[i.test]], my$projection[[i.test]])
plot(vu$projection[[i.test]], my$projection[[i.test]])
abline(0,1)

Rprof(paste0("./usr/fastFPSprof-", Sys.Date(), ".out"))
  my <- fastFPS(S, ndim=ndim, lambda=lambdas, verbose=0)
Rprof(NULL)
summaryRprof(paste0("./usr/fastFPSprof-", Sys.Date(), ".out"))


# ##############
# ## Try RcppArmadillo
# ##############
# library(inline)
# library(RcppArmadillo)
# g <- cxxfunction ( signature (vs = "numeric"),
#                    plugin = "RcppArmadillo", body = '
#                    arma::vec v = Rcpp::as < arma::vec > (vs);
#                    arma::mat op = v * v. t () ;
#                    double ip = arma::as_scalar ( v.t () * v);
#                    return Rcpp::List::create ( Rcpp::Named ("outer") =op,
#                                                   Rcpp::Named ("inner") = ip );
#                    ')
#
# g(7:11)
#
# ## Does Rcpp help?
# RreconSVD <- function(D, V){
#   return(V %*% D %*% t(V) )
# }
#
# library(microbenchmark)
# microbenchmark(R=RreconSVD(diag(D), V),
#                RcppCallfromR = reconSVDcpp(diag(D), V),
#                RcppQuote = reconSVD(diag(D), V))
# # Call from R is the slowest
#
# load_all()
# lambdas = sort(quantile(S[lower.tri(S)], probs=c(0.8, 0.9, 0.95, 0.99)), decreasing=TRUE)
# microbenchmark(Vu=fps(S, ndim=ndim, lambda=lambdas),
#                my=fastFPS(S, ndim=ndim, lambda=lambdas),
#                times=20)


