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

sparsePLSCA <- function(X, Y, Mx = NULL, My = NULL, Wx = NULL, Wy = NULL, components = 2L, tol = .Machine$double.eps,
                        doublecentering_X = TRUE, doublecentering_Y = TRUE,
                        init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                        rdsLeft = rep(1, components), rdsRight = rep(1, components),
                        grpLeft = NULL, grpRight = NULL,
                        orthogonality = "both",
                        OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                        projPriority = "orth",
                        projPriorityLeft = projPriority,
                        projPriorityRight = projPriority,
                        itermaxALS = 1000, itermaxPOCS = 1000,
                        epsALS = 1e-10, epsPOCS = 1e-10){
  if ( !is.matrix(X) ){
    X <- as.matrix(X,rownames.force = TRUE)
  }
  Y_is_missing <- missing(Y)
  if( !Y_is_missing ){
    if ( !is.matrix(Y) ){
      Y <- as.matrix(Y, rownames.force = TRUE)
    }
  }

  X_ca_preproc <- ca_preproc(X)
  Y_ca_preproc <- ca_preproc(Y)

  if (doublecentering_X) {
    X4gsvd <- X_ca_preproc$Z
  }else{
    X4gsvd <- X_ca_preproc$O
  }

  if (doublecentering_Y) {
    Y4gsvd <- Y_ca_preproc$Z
  }else{
    Y4gsvd <- Y_ca_preproc$O
  }

  if (is.null(Mx)){
    Mx <- diag((X_ca_preproc$m)^-1)
  }
  if (is.null(Wx)){
    Wx <- diag((X_ca_preproc$w)^-1)
  }
  if (is.null(My)){
    My <- diag((Y_ca_preproc$m)^-1)
  }
  if (is.null(Wy)){
    Wy <- diag((Y_ca_preproc$w)^-1)
  }

  sGSVD.res <- sparseGSVD(X = X4gsvd, Y = Y4gsvd, LW = Wx, RW = Wy, LM = Mx, RM = My,
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

  res <- spafac.out(sGSVD.res, X = X4gsvd, Y = Y4gsvd, LW = Wx, RW = Wy, LM = Mx, RM = My)
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
