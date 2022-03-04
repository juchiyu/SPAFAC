#' sparse STATIS
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param column.design
#' @param mfa.normalize TRUE (default) or FALSE.
#' @param sparse.Cmat TRUE or FALSE (default). If the C matrix (e.g., the RV matrix) is sparsified to obtain the weights.
#' @param sparse.Splus TRUE (default) or FALSE. If the compromise is sparsified.
#' @param components.Cmat the number of dimensions to extract from the C matrix (e.g., the RV matrix).
#' @param components.Splus the number of dimensions to extract from the compromise.
#' @param sparseOption "variable" or "subtable"
#' @param center For the \code{X} matrix: A parameter to pass through to \code{center} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param scale For the \code{X} matrix: A parameter to pass through to \code{scale} in \code{\link{scale}} function; either a logical value or numeric-alike vector of length equal to the number of columns of \code{X}.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#' @param init.Cmat
#' @param init.Splus
#' @param initLeft.Splus
#' @param initRight.Splus
#' @param seed
#' @param rdsLeft.Cmat
#' @param rdsRight.Cmat
#' @param rdsLeft.Splus
#' @param rdsRight.Splus
#' @param grp.Cmat
#' @param grpLeft.Splus
#' @param grpRight.Splus
#' @param orthogonality.Cmat
#' @param orthogonality.Splus
#' @param OrthSpaceLeft.Splus
#' @param OrthSpaceRight.Splus
#' @param projPriority.Cmat
#' @param projPriority.Splus
#' @param projPriorityLeft.Splus
#' @param projPriorityRight.Splus
#' @param itermaxALS.Cmat
#' @param itermaxPOCS.Cmat
#' @param itermaxALS.Splus
#' @param itermaxPOCS.Splus
#' @param epsALS.Cmat
#' @param epsPOCS.Cmat
#' @param epsALS.Splus
#' @param epsPOCS.Splus
#'
#' @return
#' @export
#'
#' @source Some lines of the function are inspired by the MExPosition package by Derek Beaton and Cherise Chin Fatt.
#' @examples

sparseSTATIS <- function(X, column.design, mfa.normalize = TRUE,
                      sparse.Cmat = FALSE, sparse.Splus = TRUE,
                      components.Cmat = 0, components.Splus = 0,
                      sparseOption = "variable",
                      center = TRUE, scale = TRUE,
                      tol = .Machine$double.eps,
                      init.Cmat = "eig", init.Splus = "svd", initLeft.Splus = NULL, initRight.Splus = NULL, seed = NULL,
                      rdsLeft.Cmat = rep(1, components.Cmat), rdsRight.Cmat = rep(1, components.Cmat),
                      rdsLeft.Splus = rep(1, components.Splus), rdsRight.Splus = rep(1, components.Splus),
                      grp.Cmat = NULL, grpLeft.Splus = NULL, grpRight.Splus = NULL,
                      orthogonality.Cmat = "loadings",
                      OrthSpaceLeft.Cmat = NULL, OrthSpaceRight.Cmat = NULL,
                      orthogonality.Splus = "loadings",
                      OrthSpaceLeft.Splus = NULL, OrthSpaceRight.Splus = NULL,
                      projPriority.Cmat = "orth",
                      projPriorityLeft.Cmat = projPriority.Cmat,
                      projPriorityRight.Cmat = projPriority.Cmat,
                      projPriority.Splus = "orth",
                      projPriorityLeft.Splus = projPriority.Splus,
                      projPriorityRight.Slus = projPriority.Splus,
                      itermaxALS.Cmat = 1000, itermaxPOCS.Cmat = 1000,
                      itermaxALS.Splus = 1000, itermaxPOCS.Splus = 1000,
                      epsALS.Cmat = 1e-10, epsPOCS.Cmat = 1e-10,
                      epsALS.Splus = 1e-10, epsPOCS.Splus = 1e-10){

  if ( !is.matrix(X) ){
    X <- as.matrix(X,rownames.force = TRUE)
  }

  if (scale == "ss1"){
    data <- projL2(scale(X, center = center, scale = FALSE))$x  # center and scale to ss = 1 before MFA-normalization
  }else{
    data <- scale(X, center = center, scale = scale) # center and scale before MFA-normalization
  }

  # Get tables into list
  count.col <- table(column.design)
  ncol.tab <- c(0, count.col)
  tab.idx <- matrix(nrow = 2, ncol = length(count.col), dimnames = list(c("from", "to"), names(count.col)))
  to_total = 0
  from_total = 1
  for (i in 1:length(count.col)) {
    from = ncol.tab[i] + from_total
    to =  ncol.tab[i+1] + to_total
    to_total = to
    from_total = from
    tab.idx[,i] <- c(from, to)
  }

  # MFA-normalization (optional)
  data.proc <- matrix(nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))
  data.proclist <- list(length(count.col))
  alpha.vec <- vector(length = length(column.design))
  if (mfa.normalize){
    tab.svd <- list()
    tab.d <- vector(length = length(count.col))
    for (i in 1:length(count.col)) {
      tab.svd[[i]] <- svd(data[, tab.idx[1,i]:tab.idx[2,i]])
      tab.d[i] <- tab.svd[[i]]$d[1] # normalize the center-and-scaled tables
      alpha.vec[tab.idx[1,i]:tab.idx[2,i]] <- 1/tab.d[i]^2
      data.proc[,tab.idx[1,i]:tab.idx[2,i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]/(tab.d[i]^2)
      data.proclist[[i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]/(tab.d[i]^2)
    }
  }else{
    alpha.vec <- rep(1, length(column.design))
    data.proc <- data
    for (i in 1:length(count.col)){
      data.proclist[[i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]
    }
  }

  # compute RV matrix an
  Cmat <- GetCmat(DATA.proc, RV = Cmat.is.RV, isCube = FALSE)
  ### masses for the RV matrix
  if (is.null(masses.Cmat)){
    masses.Cmat <- rep(1, nrow(Cmat))
  }else{
    ## check length
    if (length(masses.Cmat) != nrow(Cmat)){
      stop("The length of `masses.Cmat` does not equal the number of table (i.e., the third dimension of `DATA`).")
    }else{
      masses.Cmat <- masses.Cmat
    }
  }

  ## eigendecomposition Cmat

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

  class(sGSVD.res) <- c("sSVD", "sGSVD", "sGPCA", "MultiTab", "list")

  res <- spafac.out(sGSVD.res, X = data, LW = LW, RW = RW, tab.idx = tab.idx, column.design = column.design)
  res$X.preproc <- X
  res$alpha <- RW
  res$table.svd <- tab.svd
  res$table.d <- tab.d

  return(res)
}
