#' Group sparse Multiple Correspondence Analysis
#'
#' @param DATA the contingency table
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
#' @param components
#' @param make_data_nominal
#' @param correction4SI
#'
#' @return
#' @export
#'
#' @examples
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
