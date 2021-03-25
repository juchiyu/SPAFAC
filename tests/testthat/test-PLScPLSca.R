## test CA
library(ExPosition)
library(sGSVD)
library(GSVD)
library(GPLS)

data("wine", package = "GSVD")
plscor_results <- pls_cor(wine$objective, wine$subjective, components = 2)

spls.res <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "loadings", components = 2L, rdsLeft = rep(sqrt(ncol(wine$objective)), 2L), rdsRight = rep(sqrt(ncol(wine$subjective)), 2L))
