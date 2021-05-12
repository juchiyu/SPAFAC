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
spls.res.ls <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "both", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res.l <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "loadings", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res.s <- sparsePLSC(wine$objective, wine$subjective, orthogonality = "scores", components = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))

rownames(plscor_results$lx) <- rownames(spls.res.ls$lx)
rownames(plscor_results$ly) <- rownames(spls.res.ls$ly)

# test
test_that("sparsePLSC with no sparsifiaction gives back plain PLSC (orthogonality = both)", {
  expect_equal(abs(spls.res.ls$svd$u), abs(plscor_results$u), tolerance = tol) # both
  expect_equal(abs(spls.res.ls$svd$v), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.ls$svd$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.ls$eig), abs(plscor_results$l), tolerance = tol)
  expect_equal(abs(spls.res.ls$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.ls$fj), abs(plscor_results$fj), tolerance = tol)
  expect_equal(abs(spls.res.ls$lx), abs(plscor_results$lx), tolerance = tol)
  expect_equal(abs(spls.res.ls$ly), abs(plscor_results$ly), tolerance = tol)
})
test_that("sparsePLSC with no sparsifiaction gives back plain PLSC (orthogonality = loadings)", {
  expect_equal(abs(spls.res.l$svd$u), abs(plscor_results$u), tolerance = tol) # loadings
  expect_equal(abs(spls.res.l$svd$v), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.l$svd$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.l$eig), abs(plscor_results$l), tolerance = tol)
  expect_equal(abs(spls.res.l$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.l$fj), abs(plscor_results$fj), tolerance = tol)
  expect_equal(abs(spls.res.l$lx), abs(plscor_results$lx), tolerance = tol)
  expect_equal(abs(spls.res.l$ly), abs(plscor_results$ly), tolerance = tol)
})
test_that("sparsePLSC with no sparsifiaction gives back plain PLSC (orthogonality = scores)", {
  expect_equal(abs(spls.res.s$svd$u), abs(plscor_results$u), tolerance = tol) # scores
  expect_equal(abs(spls.res.s$svd$v), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.s$svd$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.s$eig), abs(plscor_results$l), tolerance = tol)
  expect_equal(abs(spls.res.s$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.s$fj), abs(plscor_results$fj), tolerance = tol)
  expect_equal(abs(spls.res.s$lx), abs(plscor_results$lx), tolerance = tol)
  expect_equal(abs(spls.res.s$ly), abs(plscor_results$ly), tolerance = tol)
})

# test PLSCA =================================================
data("snps.druguse", package = "GSVD")
X <- make_data_disjunctive(snps.druguse$DATA1)
Y <- make_data_disjunctive(snps.druguse$DATA2)

# get grping vector
var.grp.X <- sub("\\..*", "", colnames(X))
var.grp.Y <- sub("\\..*", "", colnames(Y))

# run PLSCA
plscacor_results <- plsca_cor(X, Y, components = 3)
# run sPLSCA
splsca.res.ls <- sparsePLSCA(X, Y, orthogonality = "both", components = 3L, rdsLeft = rep(sqrt(6), 3), rdsRight = rep(sqrt(6), 3), grpLeft = var.grp.X, grpRight = var.grp.Y)
splsca.res.l <- sparsePLSCA(X, Y, orthogonality = "loadings", components = 3L, rdsLeft = rep(sqrt(6), 3), rdsRight = rep(sqrt(6), 3), grpLeft = var.grp.X, grpRight = var.grp.Y)
splsca.res.s <- sparsePLSCA(X, Y, orthogonality = "scores", components = 3L, rdsLeft = rep(sqrt(6), 3), rdsRight = rep(sqrt(6), 3), grpLeft = var.grp.X, grpRight = var.grp.Y)

rownames(plscacor_results$lx) <- rownames(splsca.res.ls$lx)
rownames(plscacor_results$ly) <- rownames(splsca.res.ls$ly)

# test
test_that("sparsePLSCA with no sparsifiaction gives back plain PLSCA (orthogonality = both)", {
  expect_equal(abs(splsca.res.ls$svd$u), abs(plscacor_results$u), tolerance = tol)
  expect_equal(abs(splsca.res.ls$svd$v), abs(plscacor_results$v), tolerance = tol)
  expect_equal(abs(splsca.res.ls$svd$d), abs(plscacor_results$d), tolerance = tol)
  expect_equal(abs(splsca.res.ls$gsvd$p), abs(plscacor_results$p), tolerance = tol)
  expect_equal(abs(splsca.res.ls$gsvd$q), abs(plscacor_results$q), tolerance = tol)
  expect_equal(abs(splsca.res.ls$eig), abs(plscacor_results$l), tolerance = tol)
  expect_equal(abs(splsca.res.ls$fi), abs(plscacor_results$fi), tolerance = tol)
  expect_equal(abs(splsca.res.ls$fj), abs(plscacor_results$fj), tolerance = tol)
  expect_equal(abs(splsca.res.ls$lx), abs(plscacor_results$lx), tolerance = tol)
  expect_equal(abs(splsca.res.ls$ly), abs(plscacor_results$ly), tolerance = tol)
  expect_equal(length(unique(round(diag(cor(splsca.res.ls$lx, splsca.res.ls$lx.noMx)),6))), 1)
  expect_equal(diag(cor(splsca.res.ls$ly, splsca.res.ls$ly.noMy))-1, rep(0, 3), tolerance = tol)
})

test_that("sparsePLSCA with no sparsifiaction gives back plain PLSCA (orthogonality = loadings)", {
  expect_equal(abs(splsca.res.l$svd$u), abs(plscacor_results$u), tolerance = tol)
  expect_equal(abs(splsca.res.l$svd$v), abs(plscacor_results$v), tolerance = tol)
  expect_equal(abs(splsca.res.l$svd$d), abs(plscacor_results$d), tolerance = tol)
  expect_equal(abs(splsca.res.l$gsvd$p), abs(plscacor_results$p), tolerance = tol)
  expect_equal(abs(splsca.res.l$gsvd$q), abs(plscacor_results$q), tolerance = tol)
  expect_equal(abs(splsca.res.l$eig), abs(plscacor_results$l), tolerance = tol)
  expect_equal(abs(splsca.res.l$fi), abs(plscacor_results$fi), tolerance = tol)
  expect_equal(abs(splsca.res.l$fj), abs(plscacor_results$fj), tolerance = tol)
  expect_equal(abs(splsca.res.l$lx), abs(plscacor_results$lx), tolerance = tol)
  expect_equal(abs(splsca.res.l$ly), abs(plscacor_results$ly), tolerance = tol)
  expect_equal(length(unique(round(diag(cor(splsca.res.l$lx, splsca.res.l$lx.noMx)),6))), 1)
  expect_equal(diag(cor(splsca.res.l$ly, splsca.res.l$ly.noMy))-1, rep(0, 3), tolerance = tol)
})

test_that("sparsePLSCA with no sparsifiaction gives back plain PLSCA (orthogonality = scores)", {
  expect_equal(abs(splsca.res.s$svd$u), abs(plscacor_results$u), tolerance = tol)
  expect_equal(abs(splsca.res.s$svd$v), abs(plscacor_results$v), tolerance = tol)
  expect_equal(abs(splsca.res.s$svd$d), abs(plscacor_results$d), tolerance = tol)
  expect_equal(abs(splsca.res.s$gsvd$p), abs(plscacor_results$p), tolerance = tol)
  expect_equal(abs(splsca.res.s$gsvd$q), abs(plscacor_results$q), tolerance = tol)
  expect_equal(abs(splsca.res.s$eig), abs(plscacor_results$l), tolerance = tol)
  expect_equal(abs(splsca.res.s$fi), abs(plscacor_results$fi), tolerance = tol)
  expect_equal(abs(splsca.res.s$fj), abs(plscacor_results$fj), tolerance = tol)
  expect_equal(abs(splsca.res.s$lx), abs(plscacor_results$lx), tolerance = tol)
  expect_equal(abs(splsca.res.s$ly), abs(plscacor_results$ly), tolerance = tol)
  expect_equal(length(unique(round(diag(cor(splsca.res.s$lx, splsca.res.s$lx.noMx)),6))), 1)
  expect_equal(diag(cor(splsca.res.s$ly, splsca.res.s$ly.noMy))-1, rep(0, 3), tolerance = tol)
})
