## copied from DistatisR package (by Hervé Abdi)

GetCmat <- function(CubeCP, RV = TRUE){
  nI <- dim(CubeCP)[1]
  nJ <- dim(CubeCP)[3]
  CP2 <- array(CubeCP, dim = c(nI * nI, nJ))
  C <- t(CP2) %*% CP2
  if (RV) {
    laNorm <- sqrt(apply(CP2^2, 2, sum))
    C <- C/(t(t(laNorm)) %*% laNorm)
  }
  rownames(C) <- colnames(C) <- dimnames(CubeCP)[[3]]
  return(C)
}

ComputeSplus <- function(CubeCP, alpha){
  nI <- dim(CubeCP)[1]
  nJ <- dim(CubeCP)[3]
  Splus <- matrix(0, nI, nI)
  Splus <- apply(apply(CubeCP, c(1, 2), "*", t(alpha)),
                 c(2, 3), sum)
  return(Splus)
}
