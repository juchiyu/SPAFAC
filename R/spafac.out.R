spafac.out <- function(res, X, Y = NULL, LW = NULL, RW =NULL, compact = FALSE) {
  out <- list()
  if(is_sSVD(res)) {
    out$svd$d <- res$d
    out$svd$u <- res$U[,, drop = FALSE]
    out$svd$v <- res$V[,, drop = FALSE]
    out$iter <- res$iter
    out$eig <- res$d^2
    out$fi <- res$fi # t(t(res$U) * res$d)
    out$fj <- res$fj
  }

  if(is_sGSVD(res)) {
    out$gsvd$d <- res$d
    out$gsvd$p <- res$p # res$U / sqrt_LW
    out$gsvd$q <- res$q
    out$gsvd$LW <- LW
    out$gsvd$RW <- RW
    out$fi <- res$fi # t(t(LW %*% res$p) * res$d)
    out$fj <- res$fj
  }

  out$ci <- LW * (res$p)^2
  out$cj <- RW * (res$q)^2

  if(is_sPLS(res)) {
    out$lx <- X %*% res$U
    out$ly <- Y %*% res$V
    out$sx <- res$u
    out$sy <- res$v
    if(is_sGSVD(res)){
      out$lx <- res$Wx %*% out$lx
      out$ly <- res$Wy %*% out$ly
      out$sx <- res$Wx %*% out$p
      out$sy <- res$Wy %*% out$q
    }
    rownames(res$sx) <- colnames(X)
    rownames(res$sy) <- colnames(Y)
    rownames(res$lx) <- rownames(X)
    rownames(res$ly) <- rownames(Y)
  }

  return(out)
}
