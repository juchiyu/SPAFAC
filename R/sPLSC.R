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

  X <- scale(X, center = center_X, scale = scale_X)
  Y <- scale(Y, center = center_Y, scale = scale_Y)

  sSVD.res <- sparseSVD(X = X, Y = Y, k = components,
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

  ## skip using tol for now ##

  res$d <- sSVD.res$d
  res$l <- sSVD.res$d^2
  res$u <- sSVD.res$U[,,drop = FALSE]
  res$v <- sSVD.res$V[,, drop = FALSE]

  res$lx <- X %*% res$u
  res$ly <- Y %*% res$v
  res$fi <- t(t(res$u) * res$d)
  res$fj <- t(t(res$v) * res$d)
  res$iter <- sSVD.res$iter

  rownames(res$fi) <- rownames(res$u) <- colnames(X)
  rownames(res$fj) <- rownames(res$v) <- colnames(Y)

  rownames(res$lx) <- rownames(X)
  rownames(res$ly) <- rownames(Y)

  class(res) <- c("splsc", "sGSVD", "list")
  return(res)
}
