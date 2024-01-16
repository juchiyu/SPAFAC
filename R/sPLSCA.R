#' Sparse Partial Least Square Correspondence Analysis
#'
#' Performs sparse partial least square correspondence analysis (PLSCA) on two data matrices.
#'
#' @param X Data matrix with I rows and J columns.
#' @param Y Data matrix with I rows and K columns.
#' @param Mx Diagonal weight matrix for rows of X; defaults to NULL.
#' @param My Diagonal weight matrix for rows of Y; defaults to NULL.
#' @param Wx Diagonal weight matrix for columns of X; defaults to NULL.
#' @param Wy Diagonal weight matrix for columns of Y; defaults to NULL.
#' @param components The number of dimensions to extract; defaults to 2.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param doublecentering_X Logical indicating whether to apply double centering on X; defaults to TRUE.
#' @param doublecentering_Y Logical indicating whether to apply double centering on Y; defaults to TRUE.
#' @param init Initialization method; defaults to "svd".
#' @param initLeft Initial values for the left side; defaults to NULL.
#' @param initRight Initial values for the right side; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft Radii for the left side; defaults to rep(1, components).
#' @param rdsRight Radii for the right side; defaults to rep(1, components).
#' @param grpLeft Grouping vector for the left side; defaults to NULL.
#' @param grpRight Grouping vector for the right side; defaults to NULL.
#' @param orthogonality Type of orthogonality constraint; defaults to "both".
#' @param OrthSpaceLeft Orthogonal space for the left side; defaults to NULL.
#' @param OrthSpaceRight Orthogonal space for the right side; defaults to NULL.
#' @param projPriority Priority of the projection; defaults to "orth".
#' @param projPriorityLeft Priority of the left projection; defaults to projPriority.
#' @param projPriorityRight Priority of the right projection; defaults to projPriority.
#' @param itermaxALS Maximum number of ALS iterations; defaults to 1000.
#' @param itermaxPOCS Maximum number of POCS iterations; defaults to 1000.
#' @param epsALS Convergence criterion for ALS; defaults to 1e-10.
#' @param epsPOCS Convergence criterion for POCS; defaults to 1e-10.
#'
#' @return Returns an object containing the results of the sparse partial least square correspondence analysis. This object includes results from the sparse generalized singular value decomposition (sparseGSVD) applied to X and Y, and additional data specific to the analysis.
#'
#' @export
#'
#' @examples
#' # Example usage of sparsePLSCA function
#' # Assuming `X` and `Y` are data matrices
#' \dontrun{result <- sparsePLSCA(X = X, Y = Y)}
#'
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
