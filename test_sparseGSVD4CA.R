library(ExPosition)
library(tidyverse)

devtools::load_all("../sGSVD/R/sparseGSVD.R")

data(authors)
dat <- authors %>% as.matrix
N <- sum(dat)
dat <- (1/N)*dat

r <- rowSums(dat)
c <- colSums(dat)

X <- dat - r %*% t(c)

Dr <- r %>% diag
Dc <- c %>% diag

M <- Dr %>% solve
W <- Dc %>% solve

ca.res <- epCA(dat, graphs = FALSE)

spafac.res <- sparseGSVD(X = X, LW = M, RW = W, k = 2L, init = "svd",
                       rdsLeft = rep(sqrt(nrow(X)), 2), rdsRight = rep(sqrt(ncol(X)),2))

data.frame(epCA = ca.res$ExPosition.Data$fi, SPAFAC = spafac.res$fi)
data.frame(epCA = ca.res$ExPosition.Data$fj, SPAFAC = spafac.res$fj)
