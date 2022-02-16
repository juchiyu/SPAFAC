spafac.out <- function(res, X, Y = NULL, LW = NULL, RW =NULL, LM = NULL, RM = NULL, tab.idx = NULL, column.design = NULL, compact = FALSE, X4disc = NULL) {
  if ( is.matrix(LW) ){
    LW <- diag(LW)
    warning("Only the diagnonal of the LW matrix is used.")
  }
  if ( is.matrix(RW) ){
    RW <- diag(RW)
    warning("Only the diagnonal of the RW matrix is used.")
  }
  if ( is.matrix(LM) ){
    LM <- diag(LM)
    warning("Only the diagnonal of the LM matrix is used.")
  }
  if ( is.matrix(RM) ){
    RM <- diag(RM)
    warning("Only the diagnonal of the RM matrix is used.")
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
    out$ci <- res$U^2
    out$cj <- res$V^2
    out$sparsity$rdsLeft <- res$rdsLeft
    out$sparsity$rdsRight <- res$rdsRight
    out$sparsity$SI <- res$SI
  }

  if(is_sGSVD(res)) {
    out$gsvd$d <- res$d
    out$gsvd$p <- res$p # res$U / sqrt_LW
    out$gsvd$q <- res$q
    out$gsvd$LW <- LW
    out$gsvd$RW <- RW
    out$fi <- res$fi # t(t(LW %*% res$p) * res$d)
    out$fj <- res$fj
    out$ci <- LW * (res$p)^2
    out$cj <- RW * (res$q)^2
  }

  if(is_sEIGEN(res)) {
    out$eig$d <- res$d
    out$eig$l <- res$l
    out$eig$u <- res$u[,, drop = FALSE]
    out$iter <- res$iter
    out$f <- res$f # t(t(res$U) * res$d)
    out$c <- res$u^2
    out$input <- res$input
    out$sparsity$rds <- res$rds
    out$sparsity$SI <- res$SI
  }

  if(is_sGEIGEN(res)) {
    out$geig$d <- res$d
    out$geig$l <- res$l
    out$geig$p <- res$p
    out$geig$M <- LW
    out$f <- res$f # t(t(res$U) * res$d)
    out$c <- LW * (res$p)^2
    out$input <- res$input
    out$sparsity$rds <- res$rds
    out$sparsity$SI <- res$SI
  }

  if(is_sGPCA(res)) {
    out$fi <- t(t(res$p) * res$d)
    out$fj <- t(t(res$q) * res$d)
    out$ct <- apply(out$cj, 2, function(x) tapply(x, column.design, FUN = sum))
  }

  if(is_sPLS(res)) {
    out$lx <- X %*% res$U
    out$ly <- Y %*% res$V
    out$sx <- res$U
    out$sy <- res$V
    if(is_sGSVD(res)){
      out$gsvd$LM <- LM
      out$gsvd$RM <- RM
      out$lx <- (t(t(X) * LW) * sqrt(LM)) %*% res$p # sqrt_psd_matrix(LM) %*% X  %*% out$gsvd$LW %*% res$p
      out$ly <- (t(t(Y) * RW) * sqrt(RM)) %*% res$q # sqrt_psd_matrix(RM) %*% Y  %*% out$gsvd$RW %*% res$q
      # the next two are from the paper
      out$lx.noMx <- (t(t(X) * LW)) %*% res$p # X %*% out$gsvd$LW %*% res$p
      out$ly.noMy <- (t(t(Y) * RW)) %*% res$q # Y %*% out$gsvd$RW %*% res$q
      out$sx <- out$gsvd$p * LW # out$gsvd$LW %*% out$gsvd$p
      out$sy <- out$gsvd$q * RW # out$gsvd$RW %*% out$gsvd$q
      # rownames(out$lx.noMx) <- rownames(X)
      # rownames(out$ly.noMy) <- rownames(Y)
    }
    # rownames(out$sx) <- colnames(X)
    # rownames(out$sy) <- colnames(Y)
    # rownames(out$lx) <- rownames(X)
    # rownames(out$ly) <- rownames(Y)
  }

  if (is_discriminant(res)){ # so far, this only works for sDiCA
    ## get linear projector that transforms the columns
    InMat <- X4disc$X.disj.grp
    OutMat <- res$p %*% diag(res$d) %*% t(res$q)
    col.projector <- corpcor::pseudoinverse(t(InMat) %*% InMat) %*% t(InMat) %*% OutMat
    if (sum(OutMat - (InMat %*% col.projector)) > 1e-10 ){
      warning("There is no good linear approximate of the projector.")
    }

    ## profiles and row masses
    i.sup.prof <- X4disc$X.disj/rowSums(X4disc$X.disj)
    Dr_14fii <- diag(LW[X4disc$design])

    ## supplementary projections
    out$fii.unscale <- Dr_14fii %*% X4disc$X.disj %*% col.projector %*% diag(RW) %*% res$q
    n2scale <- colSums(X4disc$design.disj)[X4disc$design]
    out$fii <- out$fii.unscale * n2scale
    rownames(out$fii.unscale) <- rownames(out$fii) <- rownames(X4disc$X.disj)

    ## compute distance to the means
    dist2fi <- matrix(NA, nrow(out$fii), nrow(out$fi), dimnames = list(c(), rownames(out$fi)))
    for (i in 1:nrow(out$fi)){
      dist2fi[,i] <- sqrt(rowSums(t(t(out$fii)-out$fi[i,])^2))
    }
    min.idx <- apply(dist2fi, 1, which.min)
    assignment <- colnames(dist2fi)[min.idx]

    ## print classification results
    truth <- X4disc$design
    predict <- assignment
    out$classification$table <- table(truth, predict)
    out$classification$acc <- sum(diag(out$classification$table))/sum(out$classification$table)
    out$classification$acc.per.group <- diag(out$classification$table)/rowSums(out$classification$table)
  }

  if (is_multitab(res)){
    if (is_GSVD(res)){
      n.tab <- ncol(tab.idx)
      out$partial.fi <- array(dim = c(dim(res$fi), n.tab))
      for (i in 1:n.tab){
        out$partial.fi[,,i] <- n.tab * RW[tab.idx[1,i]] * X[,tab.idx[1,i]:tab.idx[2,i]] %*% res$q[tab.idx[1,i]:tab.idx[2,i],]
      }
    }else if (is_sEIGEN(res)){

    }
  }

  return(out)
}
