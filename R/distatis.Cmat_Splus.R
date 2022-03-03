## copied from DistatisR package (by Hervé Abdi)

GetCmat <- function(CubeCP, RV = TRUE, isCube = TRUE){
  if(isCube){
    nI <- dim(CubeCP)[1]
    nJ <- dim(CubeCP)[3]
    CP2 <- array(CubeCP, dim = c(nI * nI, nJ))
  }else{
    # if the input is a matrix
    if(is.list(CubeCP) == FALSE){
      stop("The tables need to be either in a Cube (when dimensions match) or in a list.")
    }
    GetS <- lapply(CubeCP, function(x){S <- tcrossprod(x)})
    CP2 <- simplify2array(lapply(GetS, as.vector))
  }
  C <- t(CP2) %*% CP2
  if (RV) {
    laNorm <- sqrt(apply(CP2^2, 2, sum))
    C <- C/(t(t(laNorm)) %*% laNorm)
  }
  if(isCube){
    rownames(C) <- colnames(C) <- dimnames(CubeCP)[[3]]
  }else{
    rownames(C) <- colnames(C) <- names(CubeCP)
  }
  return(C)
}

ComputeSplus <- function(CubeCP, alpha, isCube = FALSE){
  if(isCube){
    nI <- dim(CubeCP)[1]
    nJ <- dim(CubeCP)[3]
    Splus <- matrix(0, nI, nI)
    Splus <- apply(apply(CubeCP, c(1, 2), "*", t(alpha)),
                   c(2, 3), sum)
  }else{
    if(is.list(CubeCP) == FALSE){
      stop("The tables need to be either in a Cube (when dimensions match) or in a list.")
    }
    Splus.list <- mapply('*', CubeCP, alpha, SIMPLIFY = FALSE)
    Splus <- do.call("cbind", Splus.list)
  }
  return(Splus)
}
