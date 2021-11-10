#' Group sparse Discriminant Correspondence Analysis
#'
#' @param DATA the contingency table
#' @param design the design vector of the observations (i.e., the rows) of DATA
#' @param tol
#' @param make_data_nominal (TRUE or FALSE)
#' @param make_design_nominal
#' @param doublecentering
#' @param init
#' @param initLeft
#' @param initRight
#' @param seed
#' @param rdsLeft
#' @param rdsRight
#' @param grpLeft the grouping vector for the rows
#' @param grpRight the grouping vector for the columns
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
#' @param ca.type
#' @param components
#'
#' @return
#' @export
#'
#' @examples
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
    warning("ca.type is set to `mca` and you're running DiMCA.")

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
    warning("ca.type is set to `ca` and you're running DiSCA.")

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
