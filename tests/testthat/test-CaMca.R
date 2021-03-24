## test CA
library(ExPosition)
library(sGSVD)
library(GSVD)

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

gsvd.words.res <- gsvd(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW))

# run CA
ca.authors.res <- epCA(X.ca, graphs = FALSE)

sGSVD.res <- sparseGSVD(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW), k = 2L, init = "svd", rdsLeft = rep(sqrt(nrow(X.ca)), 2), rdsRight = rep(sqrt(ncol(X.ca)), 2))

sca.words.res <- sparseCA(as.matrix(X.ca), k = 2L, rdsLeft = rep(sqrt(nrow(X.ca)), 2), rdsRight = rep(sqrt(ncol(X.ca)), 2))

test_that("sparseCA gives back plain CA", {
  expect_equal(round(sca.words.res$gsvd$p,10), round(ca.authors.res$ExPosition.Data$M * ca.authors.res$ExPosition.Data$pdq$p,10))
  expect_equal(round(sca.words.res$gsvd$p,10), round(gsvd.words.res$p,10))
  expect_equal(round(sca.words.res$gsvd$q,10), round(ca.authors.res$ExPosition.Data$pdq$q,10))
  expect_equal(round(sca.words.res$gsvd$d,10), round(ca.authors.res$ExPosition.Data$pdq$Dv,10))
  expect_equal(round(sca.words.res$svd$u,10), round(gsvd.words.res$u,10))
  expect_equal(round(sca.words.res$svd$v,10), round(gsvd.words.res$v,10))
  expect_equal(round(sca.words.res$fi,10), round(ca.authors.res$ExPosition.Data$fi,10))
  expect_equal(round(sca.words.res$fj,10), round(ca.authors.res$ExPosition.Data$fj,10))
  expect_equal(round(sca.words.res$ci,10), round(ca.authors.res$ExPosition.Data$ci,10))
  expect_equal(round(sca.words.res$cj,10), round(ca.authors.res$ExPosition.Data$cj,10))
})

# test MCA ==============================
data(mca.wine)
X.mca <- mca.wine$data

mca.wine.res <- epMCA(X.mca, graphs = FALSE)

Jk <- length(sapply(X.mca, levels))
smca.words.res <- sparseMCA(as.matrix(X.mca), k = 2L, rdsLeft = rep(sqrt(nrow(X.mca)), 2), rdsRight = rep(sqrt(Jk), 2))

