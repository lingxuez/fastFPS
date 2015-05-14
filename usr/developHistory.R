library("devtools")
library("roxygen2")

# dependencies
devtools::use_package("fps")
devtools::use_package("ggplot2")

# add documentation
devtools::document()

# unit test
devtools::use_testthat()
devtools::test()

# develop
load_all()
test()

myResult = admm(S, ndim=2, lambda=0.01)
trueResult = fps(S, ndim=2, lambda=0.01)$projection[[1]]
plot(myResult, trueResult)
abline(0,1 )
