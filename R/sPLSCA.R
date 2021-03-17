#' sparse Partial Least Square Correlation
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

sparsePLSCA <- function(DATA, k = 0, tol = .Machine$double.eps,
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
                        epsALS = 1e-10, epsPOCS = 1e-10){
  sparseGSVD(X, LW, RW, k = 0, tol = .Machine$double.eps,
             init, initLeft = NULL, initRight = NULL, seed = NULL,
             rdsLeft = rep(1, R), rdsRight = rep(1, R),
             grpLeft = NULL, grpRight = NULL,
             orthogonality = "loadings",
             OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
             projPriority = "orth",
             projPriorityLeft = projPriority,
             projPriorityRight = projPriority,
             itermaxALS = 1000, itermaxPOCS = 1000,
             epsALS = 1e-10, epsPOCS = 1e-10)

}
