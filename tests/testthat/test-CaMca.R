## test CA
library(ExPosition)
library(sGSVD)
library(GSVD)

tol = 1e-15

## test CA =============================
data(authors)
X.ca <- authors

# preproc for gsvd
N <- sum(X.ca)
X.ca.proc <- 1/N * X.ca
Lv <- rowSums(X.ca.proc)
Rv <- colSums(X.ca.proc)
X.ca.proc <- X.ca.proc - Lv %*% t(Rv)
LW <- 1/Lv
RW <- 1/Rv

gsvd.authors.res <- gsvd(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW))

# run CA
ca.authors.res <- epCA(X.ca, graphs = FALSE)

sGSVD.res <- sparseGSVD(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW), k = 2L, init = "svd", rdsLeft = rep(sqrt(nrow(X.ca)), 2), rdsRight = rep(sqrt(ncol(X.ca)), 2))

sca.authors.res <- sparseCA(as.matrix(X.ca), components = 2L, rdsLeft = rep(sqrt(nrow(X.ca)), 2), rdsRight = rep(sqrt(ncol(X.ca)), 2))

test_that("sparseCA gives back plain CA", {
  expect_equal(sca.words.res$gsvd$p, ca.authors.res$ExPosition.Data$M * ca.authors.res$ExPosition.Data$pdq$p, tolerance = tol)
  expect_equal(sca.words.res$gsvd$p, gsvd.words.res$p, tolerance = tol)
  expect_equal(sca.words.res$gsvd$q, ca.authors.res$ExPosition.Data$pdq$q, tolerance = tol)
  expect_equal(sca.words.res$gsvd$d, ca.authors.res$ExPosition.Data$pdq$Dv, tolerance = tol)
  expect_equal(sca.words.res$svd$u, gsvd.words.res$u, tolerance = tol)
  expect_equal(sca.words.res$svd$v, gsvd.words.res$v, tolerance = tol)
  expect_equal(sca.words.res$fi, ca.authors.res$ExPosition.Data$fi, tolerance = tol)
  expect_equal(sca.words.res$fj, ca.authors.res$ExPosition.Data$fj, tolerance = tol)
  expect_equal(sca.words.res$ci, ca.authors.res$ExPosition.Data$ci, tolerance = tol)
  expect_equal(sca.words.res$cj, ca.authors.res$ExPosition.Data$cj, tolerance = tol)
})

# test MCA ==============================
data(mca.wine)
X.mca <- mca.wine$data

# preproc for gsvd

N.mca <- sum(tab_disjonctif(X.mca))
X.mca.proc <- 1/N.mca * tab_disjonctif(X.mca)
Lv.mca <- rowSums(X.mca.proc)
Rv.mca <- colSums(X.mca.proc)
X.mca.proc <- X.mca.proc - Lv.mca %*% t(Rv.mca)
LW.mca <- 1/Lv.mca
RW.mca <- 1/Rv.mca
rownames(X.mca.proc) <- rownames(X.mca)

mca.wine.res <- epMCA(X.mca, graphs = FALSE)

# number of levels after nominalizing
Jk <- length(unlist(sapply(X.mca, levels)))
# get grping vector
var.grp <- sub("\\..*", "", names(unlist(sapply(X.mca, levels))))

gsvd.wine.res <- gsvd(as.matrix(X.mca.proc), LW = diag(LW.mca), RW = diag(RW.mca),k = 3)

smca.wine.res <- sparseMCA(as.matrix(X.mca), components = 3L, rdsLeft = rep(sqrt(nrow(X.mca)), 3), rdsRight = rep(sqrt(Jk), 3), grpRight = var.grp)

test_that("sparseMCA gives back plain MCA", {
  expect_equal(abs(smca.wine.res$gsvd$p), abs(mca.wine.res$ExPosition.Data$M * mca.wine.res$ExPosition.Data$pdq$p), tolerance = 1e-10)
  expect_equal(abs(smca.wine.res$gsvd$p), abs(gsvd.wine.res$p), tolerance = tol)
  expect_equal(abs(smca.wine.res$gsvd$q), abs(mca.wine.res$ExPosition.Data$pdq$q[rownames(smca.wine.res$gsvd$q),]), tolerance = tol)
  expect_equal(smca.wine.res$gsvd$d, mca.wine.res$ExPosition.Data$pdq.uncor$Dv[1:3], tolerance = tol)
  expect_equal(abs(smca.wine.res$svd$u), abs(gsvd.wine.res$u), tolerance = tol)
  expect_equal(abs(smca.wine.res$svd$v), abs(gsvd.wine.res$v), tolerance = tol)
  expect_equal(abs(smca.wine.res$fi), abs(gsvd.wine.res$fi), tolerance = tol)
  expect_equal(abs(smca.wine.res$fj), abs(gsvd.wine.res$fj), tolerance = tol)
  expect_equal(abs(smca.wine.res$ci), abs(mca.wine.res$ExPosition.Data$ci), tolerance = tol)
  expect_equal(abs(smca.wine.res$cj), abs(mca.wine.res$ExPosition.Data$cj[rownames(smca.wine.res$cj),]), tolerance = tol)
})
