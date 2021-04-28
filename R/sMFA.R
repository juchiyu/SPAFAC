#' sparse Multiple Factor Analysis
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param column.design
#' @param components the number of dimensions to extract
#' @param sparseOption "variable" or "subtable
#' @param center For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#' @param init
#' @param initLeft
#' @param initRight
#' @param seed
#' @param rdsLeft
#' @param rdsRight
#' @param grpLeft
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
#' @source Some lines of the function are inspired by the MExPosition package by Derek Beaton and Cherise Chin Fatt.
#' @examples

sparseMFA <- function(X, column.design, components = 0,
                      sparseOption = "variable",
                       center = TRUE, scale = TRUE, mfa.scale = TRUE,
                       tol = .Machine$double.eps,
                       init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                       rdsLeft = rep(1, components), rdsRight = rep(1, components),
                       grpLeft = NULL,
                       orthogonality = "loadings",
                       OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                       projPriority = "orth",
                       projPriorityLeft = projPriority,
                       projPriorityRight = projPriority,
                       itermaxALS = 1000, itermaxPOCS = 1000,
                       epsALS = 1e-10, epsPOCS = 1e-10){

  if ( !is.matrix(X) ){
    X <- as.matrix(X,rownames.force = TRUE)
  }

  ### We should use the GSVD instead of the weighted PCA !!!

  data <- scale(X, center = center, scale = scale) # center and scale before MFA-normalization

  if (mfa.scale){
    alpha.vec <- vector(length = length(column.design))
    count.col <- table(column.design)
    ncol.tab <- c(0, count.col)
    tab.idx <- matrix(nrow = 2, ncol = length(count.col), dimnames = list(c("from", "to"), names(count.col)))
    tab.d <- vector(length = length(count.col))
    to_total = 0
    from_total = 1
    for (i in 1:length(count.col)) {
      from = ncol.tab[i] + from_total
      to =  ncol.tab[i+1] + to_total
      to_total = to
      from_total = from
      tab.d[i] <- svd(data[, from:to])$d[1] # normalize the center-and-scaled tables
      alpha.vec[from:to] <- 1/tab.d[i]^2
      tab.idx[,i] <- c(from, to)
    }
  }else{
    alpha.vec <- rep(1, length(column.design))
  }


  if (sparseOption == "variable"){
    grpRight = NULL
  }else if (sparseOption == "subtable"){
    grpRight = column.design
  }

  LW <- rep(1/nrow(data), nrow(data))
  RW <- alpha.vec

  sGSVD.res <- sparseGSVD(X = data, LW = LW, RW = RW, k = components,
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

  class(sGSVD.res) <- c("sSVD", "MultiTab", "list")

  res <- spafac.out(sGSVD.res, X = data, LW = LW, RW = RW, tab.idx = tab.idx)
  res$X.preproc <- X
  res$alpha <- RW
  res$table.d <- tab.d

  return(res)
}
