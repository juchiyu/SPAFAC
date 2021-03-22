#' sparse Partial Least Square Correlation
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param components the number of dimensions to extract
#' @param center_X For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param center_Y For the \code{Y} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param scale_X For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale_Y For the \code{Y} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{Y}.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
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
#' @source Some arguments of the function are inspired by the GPLS package by Derek Beaton.
#' @examples

sparsePLSC <- function(X, Y, components = 0,
                       center_X = TRUE, center_Y = TRUE,
                       scale_X = TRUE, scale_Y = TRUE,
                       tol = .Machine$double.eps,
                       init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                       rdsLeft = rep(1, components), rdsRight = rep(1, components),
                       grpLeft = NULL, grpRight = NULL,
                       orthogonality = "loadings",
                       OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                       projPriority = "orth",
                       projPriorityLeft = projPriority,
                       projPriorityRight = projPriority,
                       itermaxALS = 1000, itermaxPOCS = 1000,
                       epsALS = 1e-10, epsPOCS = 1e-10){

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

  res <- spafac.out(sSVD.res, X = X, Y = Y)
  res$X.preproc <- X4svd
  res$Y.preproc <- Y4svd

  return(res)
}
