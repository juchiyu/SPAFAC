#' Group sparse Discriminant Correspondence Analysis
#'
#' @param DATA the contingency table
#' @param design the design vector of the observations (i.e., the rows) of DATA
#' @param k the number of dimensions
#' @param tol
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
#'
#' @return
#' @export
#'
#' @examples
sparseDiCA <- function(
  DATA, design, components = 0, tol = .Machine$double.eps,
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

  if (make_data_nominal){
    DATA.disj <- tab_disjonctif(DATA)
  }else{
    DATA.disj <- DATA
  }

  if (make_design_nominal){
    design.disj <- tab_disjonctif(design)
    colnames(design.disj) <- sub("^.*\\.", "", colnames(design.disj))
  }else{
    design.disj <- design
  }
    DATA.in <- t(design.disj) %*% DATA.disj

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

  class(sGSVD.res) <- c("sSVD", "sGSVD", "list")
  res <- spafac.out(sGSVD.res, X = DATA, LW = LW, RW = RW)
  res$data$X.disj <- DATA.disj
  res$data$X.preproc <- X

  return(res)

}