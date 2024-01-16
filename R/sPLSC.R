#' Sparse Partial Least Square Correlation
#'
#' Performs sparse partial least square correlation analysis on two data matrices.
#'
#' @param X Data matrix with I rows and J columns.
#' @param Y Data matrix with I rows and K columns.
#' @param components The number of dimensions to extract; defaults to 0.
#' @param center_X Logical or numeric vector for centering each column of X; passed to \code{\link{scale}}.
#' @param center_Y Logical or numeric vector for centering each column of Y; passed to \code{\link{scale}}.
#' @param scale_X Logical or numeric vector for scaling each column of X; passed to \code{\link{scale}}.
#' @param scale_Y Logical or numeric vector for scaling each column of Y; passed to \code{\link{scale}}.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
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
#' @return Returns an object containing the results of the sparse partial least square correlation analysis. This object includes results from the sparse generalized singular value decomposition (sparseGSVD) applied to X and Y, and additional data specific to the analysis.
#'
#' @export
#'
#' @source Some arguments of the function are inspired by the GPLS package by Derek Beaton.
#' @examples
#' # Example usage of sparsePLSC function
#' # Assuming `X` and `Y` are data matrices
#' \dontrun{result <- sparsePLSC(X = X, Y = Y)}
#'
sparsePLSC <- function(X, Y, components = 0,
                       center_X = TRUE, center_Y = TRUE,
                       scale_X = TRUE, scale_Y = TRUE,
                       tol = .Machine$double.eps,
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

  X4svd <- scale(X, center = center_X, scale = scale_X)
  Y4svd <- scale(Y, center = center_Y, scale = scale_Y)

  sSVD.res <- sparseGSVD(X = X4svd, Y = Y4svd, k = components,
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

  class(sSVD.res) <- c("sPLS", "sSVD", "list")

  res <- spafac.out(sSVD.res, X = X4svd, Y = Y4svd)
  res$X.preproc <- X4svd
  res$Y.preproc <- Y4svd

  return(res)
}
