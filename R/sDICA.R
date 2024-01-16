#' Group Sparse Discriminant Correspondence Analysis
#'
#' Performs group sparse discriminant correspondence analysis (DiCA) on a given contingency table, using either multiple correspondence analysis (DiMCA) or simple correspondence analysis (DiSCA) based on the specified `ca.type`.
#'
#' @param DATA The contingency table for analysis.
#' @param design The design vector of the observations (i.e., the rows) of DATA.
#' @param ca.type Type of correspondence analysis to use, either "mca" for multiple correspondence analysis or "sca" for simple correspondence analysis; defaults to "mca".
#' @param components The number of dimensions for the analysis; defaults to 0.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param make_data_nominal Logical indicating whether to treat data as nominal (TRUE) or not (FALSE); defaults to TRUE.
#' @param make_design_nominal Logical indicating whether to treat design as nominal (TRUE) or not (FALSE); defaults to TRUE.
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
#' @return Returns an object containing the results of the group sparse discriminant correspondence analysis. The specific type of analysis and the results structure depend on the `ca.type` parameter.
#'
#' @examples
#' \dontrun{
#' res.sdisca <- sparseDiCA(DATA = FrEnAuthors.unblnc,
#'   design = data.dx.unblnc$GroupDiCA.recode,
#'   ca.type = "sca", components = kopt,
#'   make_design_nominal = TRUE,
#'   rdsLeft = rep(rdsleftopt, K),
#'   rdsRight = rep(rdsrightopt, K),
#'   initLeft = U0, initRight = V0,
#'   grpLeft = NULL, grpRight = NULL)
#' }
#'
#' @export
sparseDiCA <- function(
  DATA, design, ca.type = "mca", components = 0, tol = .Machine$double.eps,
  make_data_nominal = TRUE, make_design_nominal = TRUE,
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

  if (ca.type == "mca"){
    message(crayon::cyan("ca.type is set to `mca` and you're running DiMCA."))

    res <- sparseDiMCA(
      DATA = DATA, design = design, components = components, tol = tol,
      make_data_nominal = make_data_nominal, make_design_nominal = make_design_nominal,
      doublecentering = doublecentering,
      init = init, initLeft = initLeft, initRight = initRight, seed = seed,
      rdsLeft = rdsLeft, rdsRight = rdsRight,
      grpLeft = grpLeft, grpRight = grpRight,
      orthogonality = orthogonality,
      OrthSpaceLeft = OrthSpaceLeft, OrthSpaceRight = OrthSpaceRight,
      projPriority = projPriority,
      projPriorityLeft = projPriorityLeft,
      projPriorityRight = projPriorityRight,
      itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
      epsALS = epsALS, epsPOCS = epsPOCS)

  }else if (ca.type == "sca"){
    message(crayon::cyan("ca.type is set to `sca` and you're running DiSCA."))

    res <- sparseDiSCA(
      DATA = DATA, design = design, components = components, tol = tol,
      make_design_nominal = make_design_nominal,
      doublecentering = doublecentering,
      init = init, initLeft = initLeft, initRight = initRight, seed = seed,
      rdsLeft = rdsLeft, rdsRight = rdsRight,
      grpLeft = grpLeft, grpRight = grpRight,
      orthogonality = orthogonality,
      OrthSpaceLeft = OrthSpaceLeft, OrthSpaceRight = OrthSpaceRight,
      projPriority = projPriority,
      projPriorityLeft = projPriorityLeft,
      projPriorityRight = projPriorityRight,
      itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
      epsALS = epsALS, epsPOCS = epsPOCS)

  }else{
    stop("You need to specify ca.type by choosing between `mca`(default) and `ca`.")
  }

  return(res)

}
