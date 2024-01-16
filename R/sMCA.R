#' Group Sparse Multiple Correspondence Analysis
#'
#' Performs group sparse multiple correspondence analysis (MCA) on a given contingency table.
#'
#' @param DATA The contingency table for analysis.
#' @param components The number of dimensions for the analysis; defaults to 0.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param make_data_nominal Logical indicating whether to treat DATA as nominal (TRUE) or not (FALSE); defaults to TRUE.
#' @param doublecentering Logical indicating whether double centering should be applied; defaults to TRUE.
#' @param init Method for initialization; defaults to "svd".
#' @param initLeft Initial values for the left side; defaults to NULL.
#' @param initRight Initial values for the right side; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft Radii for the left side; defaults to rep(1, components).
#' @param rdsRight Radii for the right side; defaults to rep(1, components).
#' @param grpLeft Grouping vector for the rows; defaults to NULL.
#' @param grpRight Grouping vector for the columns; required if DATA is not treated as nominal.
#' @param orthogonality Type of orthogonality constraint; defaults to "loadings".
#' @param OrthSpaceLeft Orthogonal space for the left side; defaults to NULL.
#' @param OrthSpaceRight Orthogonal space for the right side; defaults to NULL.
#' @param projPriority Priority of the projection; defaults to "orth".
#' @param projPriorityLeft Priority of the left projection; defaults to projPriority.
#' @param projPriorityRight Priority of the right projection; defaults to projPriority.
#' @param correction4SI Correction method for supplementary items; defaults to "mca".
#' @param itermaxALS Maximum number of ALS iterations; defaults to 1000.
#' @param itermaxPOCS Maximum number of POCS iterations; defaults to 1000.
#' @param epsALS Convergence criterion for ALS; defaults to 1e-10.
#' @param epsPOCS Convergence criterion for POCS; defaults to 1e-10.
#'
#' @return Returns an object of class c("sSVD", "sGSVD", "list") containing the results of the sparse multiple correspondence analysis. This object includes the results of sparse generalized singular value decomposition (sparseGSVD) and additional data specific to MCA.
#'
#' @export
#'
#' @examples
#' # Example usage of sparseMCA function
#' # Assuming `data` is a contingency table
#' #' \dontrun{result <- sparseMCA(DATA = data)}
sparseMCA <- function(
  DATA, components = 0, tol = .Machine$double.eps,
  make_data_nominal = TRUE, doublecentering = TRUE,
  init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
  rdsLeft = rep(1, components), rdsRight = rep(1, components),
  grpLeft = NULL, grpRight = NULL,
  orthogonality = "loadings",
  OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
  projPriority = "orth",
  projPriorityLeft = projPriority,
  projPriorityRight = projPriority,
  correction4SI = "mca",
  itermaxALS = 1000, itermaxPOCS = 1000,
  epsALS = 1e-10, epsPOCS = 1e-10) {

  # add make_data_nominal
  if (make_data_nominal){
    DATA.disj <- tab_disjonctif(as.matrix(DATA))
    rownames(DATA.disj) <- rownames(DATA)
    if (is.null(grpRight)){
      # create grpRight
      grpRight <- sub("\\..+","", colnames(DATA.disj))
    }else{
      # check that group size is in the correct format
      if (length(grpRight) != ncol(DATA.disj)) stop("The length of grpRight doesn't match the data. Check your grpRight!")
    }
  }else{
    DATA.disj <- as.matrix(DATA)
    if (is.null(grpRight)){
      stop("Please specify grpRight when the input data is in disjunctive code.")
    }else{
      if (length(grpRight) != ncol(DATA.disj)) stop("The length of grpRight doesn't match the data. Check your grpRight!")
    }
  }

  N <- sum(DATA.disj)
  X <- 1/N * DATA.disj
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
                          correction4SI = correction4SI,
                          itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
                          epsALS = epsALS, epsPOCS = epsPOCS)

  class(sGSVD.res) <- c("sSVD", "sGSVD", "list")
  res <- spafac.out(sGSVD.res, X = DATA, LW = LW, RW = RW)
  # ADD CORRECTION FOR EIGENVALUES
  res$data$X.disj <- DATA.disj
  res$data$X.preproc <- X

  return(res)

}
