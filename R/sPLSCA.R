#' sparse Partial Least Square Correspondence anlaysis
#'
#' @param DATA the contingency table
#' @param k the number of dimensions
#' @param doublecentering_X
#' @param doublecentering_Y
#' @param tol
#' @param init
#' @param initLeft
#' @param initRight
#' @param seed
#' @param rdsLeft
#' @param rdsRight
#' @param grpLeft
#' @param grpRight
#' @param orthogonality
#' @param OrthSpaceLeft
#' @param OrthSpaceRight
#' @param projPriority
#' @param projPriorityLeft
#' @param projPriorityRight
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#'
#' @return
#' @export
#'
#' @examples

sparsePLSCA <- function(X, Y, k = 0, tol = .Machine$double.eps,
                        doublecentering_X = TRUE, doublecentering_Y = TRUE,
                        init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                        rdsLeft = rep(1, k), rdsRight = rep(1, k),
                        grpLeft = NULL, grpRight = NULL,
                        orthogonality = "loadings",
                        OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                        projPriority = "orth",
                        projPriorityLeft = projPriority,
                        projPriorityRight = projPriority,
                        itermaxALS = 1000, itermaxPOCS = 1000,
                        epsALS = 1e-10, epsPOCS = 1e-10){

  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)

  X4gsvd <- X_ca_preproc$Z
  Y4gsvd <- Y_ca_preproc$Z

  Mx <- diag((X_ca_preproc$m)^-1)
  Wx <- diag((X_ca_preproc$w)^-1)
  My <- diag((Y_ca_preproc$m)^-1)
  Wy <- diag((Y_ca_preproc$w)^-1)

  sGSVD.res <- sparseGSVD(X = X4gsvd, Y = Y4gsvd, LW = Wx, RW = Wy,
             k = components,
             init = init, initLeft = initLeft, initRight = initRight, seed = seed,
             rdsLeft = rdsLeft, rdsRight = rdsRight,
             grpLeft = grpLeft, grpRight = grpRight,
             orthogonality = orthogonality,
             OrthSpaceLeft = OrthSpaceLeft, OrthSpaceRight = OrthSpaceRight,
             projPriority = projPriority,
             projPriorityLeft = projPriority,
             projPriorityRight = projPriority,
             itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
             epsALS = epsALS, epsPOCS = epsPOCS)


  ## skip using tol for now ##

  # 1) create a function for the output
  # 2) add a "compact" argument for the output function : keep p, d, q, Wx, Mx, Wy, My

  res$d <- sGSVD.res$d
  res$l <- sGSVD.res$d^2
  res$u <- sGSVD.res$U[,,drop = FALSE]
  res$v <- sGSVD.res$V[,, drop = FALSE]
  res$pdq$p <- sGSVD.res$p # put this in a pdq list
  res$pdq$q <- sGSVD.res$q

  res$lx <- X4gsvd %*% Wx %*% res$u
  res$ly <- Y4gsvd %*% Wy %*% res$v
  res$sx <- Wx %*% res$p
  res$sy <- Wy %*% res$q
  res$fi <- Wx %*% res$p %*% res$d
  res$fj <- Wy %*% res$q %*% res$d
  res$iter <- sGSVD.res$iter

  rownames(res$fi) <- rownames(res$u) <- colnames(X)
  rownames(res$fj) <- rownames(res$v) <- colnames(Y)

  rownames(res$lx) <- rownames(X)
  rownames(res$ly) <- rownames(Y)

  d, u, v, p, q, Wx, Wy, X, Y

  class(res) <- c("sPLS", "sSVD", "sGSVD", "list")
  return(res)

}

ca_preproc <- function(DATA){ # Derek's function in GPLS
  O <- DATA/sum(DATA)
  m <- rowSums(O)
  w <- colSums(O)
  E <- m %o% w
  Z <- O - E
  return(list(Z = Z, O = O, E = E, m = m, w = w))
}
