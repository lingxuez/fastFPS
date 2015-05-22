library("devtools")
library("roxygen2")
library(fps)
library(microbenchmark)
devtools::install_github("hadley/lineprof")

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
p = 150
set.seed(200)
X = matrix(rnorm(p*10), nrow=10, ncol=p)
S = cov(X)
ndim = 2
lambdas = quantile(S[lower.tri(S)], probs=c(0.99, 0.95, 0.9, 0.85, 0.8))

checkFpsEqual <- function(values){
  all.equal(values[[1]]$projection, values[[2]]$projection,
            tolerance=1e-2)
}

load_all()
test()
microbenchmark(vu <- fps(S, ndim=ndim, lambda=lambdas, maxiter=500),
               my <- fastFPS(S, ndim=ndim, lambda=lambdas, maxiter=500),
               my.lazy <- fastFPS.lazyscreen(S, ndim=ndim, lambda=lambdas, maxiter=500),
               times=2)
microbenchmark(vu <- fps(S, ndim=ndim, lambda=lambdas[3], verbose=0),
               my <- fastFPS(S, ndim=ndim, lambda=lambdas[3], verbose=0),
               check=checkFpsEqual, times=2)


microbenchmark(x <- findActivePrep(S, ndim),
               for (lambda in lambdas){x <- findActive(S, ndim, lambda)}
  )

i.test = 5
all.equal(vu$projection[[i.test]], my$projection[[i.test]])
plot(vu$projection[[i.test]], my$projection[[i.test]])
abline(0,1)

# vuSeq <- fps(S, ndim=ndim, lambda=lambdas, verbose=3, maxiter=500)
# vuSingle <- fps(S, ndim=ndim, lambda=lambdas[5], verbose=3, maxiter=500)
# all.equal(vuSeq$projection[[5]], vuSingle$projection[[1]])


##############
# ## Try RcppArmadillo
# ##############
load_all()
Rprof(paste0("./usr/fastFPSprof-", Sys.Date(), ".out"))
  my <- fastFPS(S, ndim=ndim, lambda=lambdas, verbose=0)
Rprof(NULL)
summaryRprof(paste0("./usr/fastFPSprof-", Sys.Date(), ".out"))

library(lineprof)
l <- lineprof(fastFPS(S, ndim=ndim, lambda=lambdas, verbose=0))
l

microbenchmark(eigen(S, symmetric=TRUE, only.values=TRUE),
               eigs(S, 10, symmetric=TRUE))

system.time({
  eigen(S, symmetric=TRUE, only.values=TRUE)
})
# 0.7 seconds for p=1000

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

load_all()
myInit0 <- fastFPS(S, ndim=ndim, lambda=lambdas[5], maxiter=500, verbose=3)
load_all()
myInit1 <- fastFPS(S, ndim=ndim, lambda=lambdas[5], maxiter=1000, verbose=3)
all.equal(myInit0$projection[[1]], myInit1$projection[[1]])
plot(myInit0$projection[[1]], myInit1$projection[[1]])

system.time({
  vu <- fps(S, ndim=ndim, lambda=lambdas, maxiter=500)
})

system.time({
  my <- fastFPS(S, ndim=ndim, lambda=lambdas, maxiter=500)
})

system.time({
  my.lazy <- fastFPS.lazyscreen(S, ndim=ndim, lambda=lambdas, maxiter=500)
})

system.time({
  my.cpp <- fastFPS_cpp(S, ndim=ndim, lambda=lambdas, maxiter=500)
})

all.equal(my$projection, my.lazy$projection, tolerance=0.01)
all.equal(my$projection, my.cpp$projection, tolerance=0.01)
all.equal(my.lazy$projection, my.cpp$projection)
all.equal(vu$projection, my.lazy$projection, 0.001)


##############
## Try Rcpp
##############
library(Rcpp)
cppFunction("int add(int x, int y, int z){
  int sum = x + y + z;
  return sum;
            }")

add
add(1,2,3)

lambda = lambdas[3]
microbenchmark(soft <- softThreshold(S, lambda),
               softcpp <- softThresholdCpp(S, lambda))
all.equal(soft, softcpp)

## workflow to build Rcpp functions:
## cmd+shift+D: document
## build&load
## try it
