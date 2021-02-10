#' Group sparse Correspondence Analysis
#'
#' @param DATA the contingency table
#' @param k the number of dimensions
#' @param tol
#' @param doublecentering
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
sparseCA <- function(
              DATA, k = 0, tol = .Machine$double.eps,
              doublecentering = TRUE,
              init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
              rdsLeft = rep(1, R), rdsRight = rep(1, R),
              grpLeft = NULL, grpRight = NULL,
              orthogonality = "loadings",
              OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
              projPriority = "orth",
              projPriorityLeft = projPriority,
              projPriorityRight = projPriority,
              itermaxALS = 1000, itermaxPOCS = 1000,
              epsALS = 1e-10, epsPOCS = 1e-10) {


  mRP <- ExPosition::makeRowProfiles(DATA, weights = NULL, masses = NULL, hellinger = FALSE)
  # N <- sum(DATA)
  # X <- 1/N * DATA
  # Lv <- rowSums(X)
  # Rv <- colSums(X)
  # if (doublecentering) X <- X - Lv %*% t(Rv)
  # LW <- 1/sqrt(Lv)
  # RW <- 1/sqrt(Rv)

  res.spgsvd <- sparseGSVD(mRP$deviations, LW = mRP$masses, RW = mRP$weights, k = k, tol = .Machine$double.eps,
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


}
