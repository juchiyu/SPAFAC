library(ExPosition)
library(sGSVD)
library(GSVD)
library(GPLS)

tol = 1e-14

# test sPLSC =============================

data("wine", package = "GSVD")

# run PLSC
plscor_results <- pls_cor(wine$objective, wine$subjective, components = 0)
# run sPLSC with no sparsification
spls.res <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "both", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "loadings", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "score", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))

rownames(plscor_results$lx) <- rownames(spls.res$lx)
rownames(plscor_results$ly) <- rownames(spls.res$ly)

# test
test_that("sparsePLSC with no sparsifiaction gives back plain PLSC", {
  expect_equal(abs(spls.res$svd$u), abs(plscor_results$u), tolerance = tol)
  expect_equal(abs(spls.res$svd$v), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res$svd$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res$svd$u), abs(plscor_results$u), tolerance = tol)
  expect_equal(abs(spls.res$eig), abs(plscor_results$l), tolerance = tol)
  expect_equal(abs(spls.res$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res$fj), abs(plscor_results$fj), tolerance = tol)
  expect_equal(abs(spls.res$lx), abs(plscor_results$lx), tolerance = tol)
  expect_equal(abs(spls.res$ly), abs(plscor_results$ly), tolerance = tol)
})

# test PLSCA =================================================
data("snps.druguse", package = "GSVD")
X <- make_data_disjunctive(snps.druguse$DATA1)
Y <- make_data_disjunctive(snps.druguse$DATA2)

# run PLSCA
plscacor_results <- plsca_cor(X, Y)
# run sPLSCA
splsca.res <- sparsePLSCA(X, Y, orthogonality = "both", components = 5L, rdsLeft = rep(sqrt(6), 5), rdsRight = rep(sqrt(6), 5))
