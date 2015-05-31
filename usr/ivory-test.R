## Test/compare on ivory module
##

## generate data
# WGCNA.controls <- read.table(paste0("~/Documents/Thesis/Data/WGCNA/",
#                                     "DLPFC.ensembl.coexpr.Control.clustering_modules.tsv"),
#                              header=TRUE)
# cpm.res <- read.csv("~/Documents/Thesis/Data/CommonMinds/CMC_SCZ_CONTROLS_RESIDUALS_16423_ENSMBL_genes_537_samples_without_SV_03-25-2015.csv",
#                     header=TRUE)
# i.cases = grep("_S", colnames(cpm.res))
# i.controls = grep("_C", colnames(cpm.res))
#
# module = "ivory"
# i.module = which(WGCNA.controls[,"Module"] == module)
# X.case = t(cpm.res[i.module, i.cases])
# colnames(X.case) = cpm.res[i.module, "GeneFeature"]
# X.control = t(cpm.res[i.module, i.controls])
# colnames(X.control) = cpm.res[i.module, "GeneFeature"]
# D.cor = cor(X.case) - cor(X.control)
#
# save(D.cor, X.case, X.control, module,
#      file=paste0("./data/", module,"Corr.RData"))

rm(list=ls())
data(ivoryCorr)

load_all()
lambdas = quantile(abs(D.cor)[lower.tri(D.cor)], probs=c(0.99, 0.95, 0.9, 0.85, 0.8))
ndim=1
system.time({
  vuFPS <- fps(D.cor, ndim=ndim, lambda=lambdas)
})
system.time({
  myFPS.lazy <- fastFPS.lazyscreen(D.cor, ndim=ndim, lambda=lambdas)
})
system.time({
  myFPS.cpp <- fastFPS.cpp(D.cor, ndim=ndim, lambda=lambdas)
})
system.time({
  myFPS <- fastFPS(D.cor, ndim=ndim, lambda=lambdas)
})

plot(vuFPS$projection[[3]], myFPS$projection[[3]])
abline(0,1)
sort(vuFPS$projection[[1]], decreasing=TRUE)[1:5]
plot(vuFPS$leverage[,1], myFPS.lazy$leverage[,1])
abline(0,1)


## Is it because D is not positive definite?
Deig <- eigen(D)
i.pos = which(Deig$values > 0.01)
Dpos <- Deig$vectors[, i.pos] %*% diag(Deig$values[i.pos]) %*% t(Deig$vectors[, i.pos])
plot(eigen(Dpos)$values)

system.time({
  vuFPS <- fps(Dpos, ndim=ndim, lambda=lambdas[2:5])
})
system.time({
  myFPS.lazy <- fastFPS.lazyscreen(Dpos, ndim=ndim, lambda=lambdas[2:5])
})
system.time({
  myFPS <- fastFPS(Dpos, ndim=ndim, lambda=lambdas[2:5])
})
system.time({
  myFPS.cpp <- fastFPS_cpp(Dpos, ndim=ndim, lambda=lambdas[2:5])
})

plot(vuFPS$projection[[4]], myFPS.cpp$projection[[4]])
abline(0,1)
all.equal(vuFPS$projection, myFPS$projection)


system.time({
  vuFPS <- fps(cor(X.case), ndim=ndim, lambda=lambdas[2:5], verbose=3)
})
system.time({
  myFPS.lazy <- fastFPS.lazyscreen(cor(X.case), ndim=ndim, lambda=lambdas[2:5], verbose=3)
})
system.time({
  myFPS.cpp <- fastFPS.cpp(cor(X.case), ndim=ndim, lambda=lambdas[2:5])
})
plot(vuFPS$projection[[1]], myFPS.cpp$projection[[1]])

v = sort(rnorm(100, sd=0.1), decreasing=TRUE)
microbenchmark(mytheta <- getTheta(v, 20),
               vusimplex <- simplex(v, 20))
