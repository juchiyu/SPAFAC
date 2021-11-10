# library(MExPosition)
library(FactoMineR)
library(sGSVD)
library(SPAFAC)

tol = 1e-14

## plain MFA ----
data("wines2007")
facto.res <- MFA(wines2007$data, group=table(as.matrix(wines2007$table)),
           ncp=5, name.group=c("E1","E2","E3"), graph = FALSE)

# mexpo.res <- mpMFA(wines2007$data, wines2007$table, graphs = FALSE)

## sparse MFA ----
data.mfa <- wines2007$data
col.design <- wines2007$table

k = 5L
rdsLeft = rep(sqrt(nrow(data.mfa)), k)
rdsRight = rep(sqrt(ncol(data.mfa)), k)

smfa.res <- sparseMFA(X = as.matrix(data.mfa), column.design = as.matrix(col.design),
                      components = k,
                      sparseOption = "variable",
                      center = TRUE, scale = TRUE, mfa.scale = TRUE,
                      orthogonality = "loadings",
                      rdsLeft = rdsLeft, rdsRight = rdsRight,
                      grpLeft = NULL)

## check results ----
facto.res$eig[,1]/smfa.res$eig # off by # of rows
# mexpo.res$mexPosition.Data$Table$eigs/smfa.res$eig # off by 1/(# of tables)
facto.res$ind$coord/smfa.res$fi # off by 1/sqrt(LW)
# mexpo.res$mexPosition.Data$Table$fi/smfa.res$fi # off by 1/sqrt(# of tables)
facto.res$ind$coord/mexpo.res$mexPosition.Data$Table$fi # off by sqrt(# of rows * # of tables)
apply(smfa.res$partial.fi, c(1,2), mean)/smfa.res$fi

facto.res$quanti.var$coord/smfa.res$fj # off by an unknown scaling factor
# mexpo.res$mexPosition.Data$Table$Q/smfa.res$gsvd$q # off by sqrt(# of tables/RW)

## test ----
test_that("sparseMFA with no sparsification gives back plain MFA", {
  expect_equivalent(abs(facto.res$eig[,1]), abs(smfa.res$eig*nrow(data.mfa)), tolerance = tol)
  expect_equivalent(abs(facto.res$ind$coord), abs(smfa.res$fi*(1/sqrt(smfa.res$gsvd$LW))), tolerance = tol)
  expect_equivalent(apply(smfa.res$partial.fi, c(1,2), mean), smfa.res$fi)
})
