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

  X4gsvd <- ca_preproc(X)$Z
  Y4gsvd <- ca_preproc(Y)$Z

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
  class(sGSVD.res) <- c("sPLS", "sSVD", "sGSVD", "list")

  res <- spafac.out(sGSVD.res, X = X, Y = Y, LW = Wx, RW = Wy)
  res$X.preproc <- X4gsvd
  res$Y.preproc <- Y4gsvd

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
