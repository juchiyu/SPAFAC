spafac.out <- function(res, X, Y = NULL, LW = NULL, RW =NULL, LM = NULL, RM = NULL, compact = FALSE) {
  if ( is.vector(LW) ){
    LW <- diag(LW)
  }
  if ( is.vector(RW) ){
    RW <- diag(RW)
  }
  out <- list()
  if(is_sSVD(res)) {
    out$svd$d <- res$d
    out$svd$u <- res$U[,, drop = FALSE]
    out$svd$v <- res$V[,, drop = FALSE]
    out$iter <- res$iter
    out$eig <- res$d^2
    out$fi <- res$fi # t(t(res$U) * res$d)
    out$fj <- res$fj
    out$ci <- res$fi^2
    out$cj <- res$fj^2
  }

  if(is_sGSVD(res)) {
    out$gsvd$d <- res$d
    out$gsvd$p <- res$p # res$U / sqrt_LW
    out$gsvd$q <- res$q
    out$gsvd$LW <- LW
    out$gsvd$RW <- RW
    out$fi <- res$fi # t(t(LW %*% res$p) * res$d)
    out$fj <- res$fj
    out$ci <- LW %*% (res$p)^2
    out$cj <- RW %*% (res$q)^2
  }



  if(is_sPLS(res)) {
    out$lx <- X %*% res$U
    out$ly <- Y %*% res$V
    out$sx <- res$U
    out$sy <- res$V
    if(is_sGSVD(res)){
      out$gsvd$LM <- LM
      out$gsvd$RM <- RM
      out$lx <- sqrt_psd_matrix(LM) %*% X  %*% out$gsvd$LW %*% res$p
      out$ly <- sqrt_psd_matrix(RM) %*% Y  %*% out$gsvd$RW %*% res$q
      out$lx.noMx <- X %*% out$gsvd$LW %*% res$p # this is from the paper
      out$ly.noMy <- Y %*% out$gsvd$RW %*% res$q
      out$sx <- out$gsvd$LW %*% out$gsvd$p
      out$sy <- out$gsvd$RW %*% out$gsvd$q
    }
    rownames(out$sx) <- colnames(X)
    rownames(out$sy) <- colnames(Y)
    rownames(out$lx) <- rownames(X)
    rownames(out$ly) <- rownames(Y)
    rownames(out$lx.noMx) <- rownames(X)
    rownames(out$ly.noMy) <- rownames(Y)
  }

  return(out)
}
