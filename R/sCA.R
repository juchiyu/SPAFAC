#' Group sparse Correspondence Analysis
#'
#' @param DATA The contingency table for analysis.
#' @param components The number of dimensions for the analysis; defaults to 0.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param doublecentering Logical indicating whether double centering should be applied; defaults to TRUE.
#' @param init Method for initialization; defaults to "svd".
#' @param initLeft Initial values for the left side; defaults to NULL.
#' @param initRight Initial values for the right side; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft Radii for the left side; defaults to rep(1, components).
#' @param rdsRight Radii for the right side; defaults to rep(1, components).
#' @param grpLeft Group assignments for the left side; defaults to NULL.
#' @param grpRight Group assignments for the right side; defaults to NULL.
#' @param orthogonality Type of orthogonality constraint; defaults to "loadings".
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
#'
#' @return Returns an object of class c("sSVD", "sGSVD", "list") containing the results
#' of the group sparse correspondence analysis. This object includes various components
#' produced by the sparse generalized singular value decomposition (sparseGSVD) and
#' the subsequent `spafac.out` processing. Key components include the singular values,
#' left and right singular vectors, and additional metrics relevant to the analysis.
#' The object also contains `X.preproc`, the preprocessed data matrix used in the analysis.
#'
#' @examples
#'
#' data("example1_sCA") # Load the Cause of Deaths Data
#' sca.res <- sparseCA(as.matrix(example1_sCA), components = 2L, rdsLeft = c(1.3, 1.3), rdsRight = c(2.3, 2.3))
#'
#' @export
sparseCA <- function(
  DATA, components = 0, tol = .Machine$double.eps,
  doublecentering = TRUE,
  init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
  rdsLeft = rep(1, components), rdsRight = rep(1, components),
  grpLeft = NULL, grpRight = NULL,
  orthogonality = "loadings",
  OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
  projPriority = "orth",
  projPriorityLeft = projPriority,
  projPriorityRight = projPriority,
  itermaxALS = 1000, itermaxPOCS = 1000,
  epsALS = 1e-10, epsPOCS = 1e-10) {


  # mRP <- ExPosition::makeRowProfiles(DATA, weights = NULL, masses = NULL, hellinger = FALSE)
  N <- sum(DATA)
  X <- 1/N * DATA
  Lv <- rowSums(X)
  Rv <- colSums(X)
  if (doublecentering) X <- X - Lv %*% t(Rv)
  LW <- 1/Lv
  RW <- 1/Rv

  sGSVD.res <- sparseGSVD(X, LW = LW, RW = RW, k = components, tol = .Machine$double.eps,
                          init = init, initLeft = initLeft, initRight = initRight, seed = NULL,
                          rdsLeft = rdsLeft, rdsRight = rdsRight,
                          grpLeft = grpLeft, grpRight = grpRight,
                          orthogonality = orthogonality,
                          OrthSpaceLeft = OrthSpaceLeft,
                          OrthSpaceRight = OrthSpaceRight,
                          projPriority = projPriority,
                          projPriorityLeft = projPriorityLeft,
                          projPriorityRight = projPriorityRight,
                          itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
                          epsALS = epsALS, epsPOCS = epsPOCS)

  class(sGSVD.res) <- c("sSVD", "sGSVD", "list")
  res <- spafac.out(sGSVD.res, X = DATA, LW = LW, RW = RW)
  res$X.preproc <- X

  return(res)

}
