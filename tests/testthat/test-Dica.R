library(TExPosition)
library(sGSVD)
library(GSVD)

tol = 1e-14

## test CA =============================
data(dica.wine)

X.data <- as.matrix(dica.wine$data)
X.des <- as.matrix(dica.wine$design)

X.dica <- t(X.des) %*% X.data

# preproc for gsvd
N <- sum(X.dica)
X.dica.proc <- 1/N * X.dica
Lv <- rowSums(X.dica.proc)
Rv <- colSums(X.dica.proc)
X.dica.proc <- X.dica.proc - Lv %*% t(Rv)
LW <- 1/Lv
RW <- 1/Rv

grpRight.constraint <- sub("\\..+","", colnames(X.data))

# run DiCA
dica.res <- tepDICA(dica.wine$data,DESIGN=dica.wine$design,make_design_nominal=FALSE, graphs = FALSE)
# run gsvd
gsvd.wine.res <- gsvd(as.matrix(X.dica.proc), LW = diag(LW), RW = diag(RW))
# run sGSVD
sGSVD.res <- sparseGSVD(as.matrix(X.dica.proc), LW = diag(LW), RW = diag(RW), k = 2L, init = "svd", rdsLeft = rep(sqrt(nrow(X.dica)), 2), rdsRight = rep(sqrt(ncol(X.dica)), 2))
# run sDiCA
sdica.wine.res <- sparseDiCA(as.matrix(X.data), design = X.des, make_data_nominal = FALSE, make_design_nominal = FALSE, components = 2L, rdsLeft = rep(sqrt(nrow(X.dica)), 2), rdsRight = rep(sqrt(ncol(X.dica)), 2), grpRight = grpRight.constraint)

test_that("sparseDiCA with no sparsification gives back plain DiCA", {
  expect_equal(abs(sdica.wine.res$gsvd$p), abs(dica.res$TExPosition.Data$M * dica.res$TExPosition.Data$pdq$p), tolerance = tol)
  expect_equal(abs(sdica.wine.res$gsvd$p), abs(gsvd.wine.res$p), tolerance = tol)
  expect_equal(abs(sdica.wine.res$gsvd$q), abs(dica.res$TExPosition.Data$pdq$q), tolerance = tol)
  expect_equal(abs(sdica.wine.res$gsvd$d), abs(dica.res$TExPosition.Data$pdq$Dv), tolerance = tol)
  expect_equal(abs(sdica.wine.res$svd$u), abs(gsvd.wine.res$u), tolerance = tol)
  expect_equal(abs(sdica.wine.res$svd$v), abs(gsvd.wine.res$v), tolerance = tol)
  expect_equal(abs(sdica.wine.res$fi), abs(dica.res$TExPosition.Data$fi), tolerance = tol)
  expect_equal(abs(sdica.wine.res$fj), abs(dica.res$TExPosition.Data$fj), tolerance = tol)
  expect_equal(abs(sdica.wine.res$ci), abs(dica.res$TExPosition.Data$ci), tolerance = tol)
  expect_equal(abs(sdica.wine.res$cj), abs(dica.res$TExPosition.Data$cj), tolerance = tol)
})
