#' Group Sparse Discriminant Correspondence Analysis for Simple CA
#'
#' Performs group sparse discriminant correspondence analysis specifically for simple correspondence analysis (SCA) on a given contingency table.
#'
#' @param DATA The contingency table for analysis.
#' @param design The design vector of the observations (rows) of DATA.
#' @param components The number of dimensions for the analysis; defaults to 0.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param make_design_nominal Logical indicating whether to treat design vector as nominal (TRUE) or not (FALSE); defaults to TRUE.
#' @param doublecentering Logical indicating whether double centering should be applied; defaults to TRUE.
#' @param init Method for initialization; defaults to "svd".
#' @param initLeft Initial values for the left side; defaults to NULL.
#' @param initRight Initial values for the right side; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft Radii for the left side; defaults to rep(1, components).
#' @param rdsRight Radii for the right side; defaults to rep(1, components).
#' @param grpLeft Grouping vector for the rows; defaults to NULL.
#' @param grpRight Grouping vector for the columns; defaults to NULL.
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
#' @return Returns an object of class c("sSVD", "sGSVD", "SCA", "discriminant", "list") containing the results of the analysis. This object includes the results of sparse generalized singular value decomposition (sparseGSVD) and additional data specific to discriminant SCA.
#'
#' @export
#'
#' @examples
#' # Example usage of sparseDiSCA function
#' # Assuming `data` is a contingency table and `design` is the design vector
#' \dontrun{result <- sparseDiSCA(DATA = data, design = design)}
#'
sparseDiSCA <- function(
  DATA, design, components = 0, tol = .Machine$double.eps,
  make_design_nominal = TRUE,
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

  if (make_design_nominal){
    design.disj <- tab_disjonctif(design)
    colnames(design.disj) <- sub("^.*\\.", "", colnames(design.disj))
  }else{
    design.disj <- design
    design <- colnames(design.disj)[apply(design.disj, 1, function(x)which(x>0))]
  }

  DATA.in <- t(design.disj) %*% DATA

  N <- sum(DATA.in)
  X <- 1/N * DATA.in
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

  class(sGSVD.res) <- c("sSVD", "sGSVD", "SCA", "discriminant", "list")
  DATA.disj <- DATA
  res <- spafac.out(sGSVD.res, X = X, LW = LW, RW = RW, X4disc = list(design = design, design.disj = design.disj, X.disj = DATA.disj, X.disj.grp = DATA.in))
  res$data$X.disj <- DATA.disj
  res$data$X.preproc <- X

  return(res)

}
