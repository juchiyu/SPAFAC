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
                      sparse.Cmat = FALSE, sparse.grandX = TRUE,
                      components.Cmat = 0, components.grandX = 0,
                      sparseOption = "variable",
                      center = TRUE, scale = TRUE,
                      tol = .Machine$double.eps,
                      init.Cmat = "eig", init.grandX = "svd", initLeft.grandX = NULL, initRight.grandX = NULL, seed = NULL,
                      rdsLeft.Cmat = rep(1, components.Cmat), rdsRight.Cmat = rep(1, components.Cmat),
                      rdsLeft.grandX = rep(1, components.grandX), rdsRight.grandX = rep(1, components.grandX),
                      grp.Cmat = NULL, grpLeft.grandX = NULL, grpRight.grandX = NULL,
                      orthogonality.Cmat = "loadings",
                      OrthSpaceLeft.Cmat = NULL, OrthSpaceRight.Cmat = NULL,
                      orthogonality.grandX = "loadings",
                      OrthSpaceLeft.grandX = NULL, OrthSpaceRight.grandX = NULL,
                      projPriority.Cmat = "orth",
                      projPriorityLeft.Cmat = projPriority.Cmat,
                      projPriorityRight.Cmat = projPriority.Cmat,
                      projPriority.grandX = "orth",
                      projPriorityLeft.grandX = projPriority.grandX,
                      projPriorityRight.Slus = projPriority.grandX,
                      itermaxALS.Cmat = 1000, itermaxPOCS.Cmat = 1000,
                      itermaxALS.grandX = 1000, itermaxPOCS.grandX = 1000,
                      epsALS.Cmat = 1e-10, epsPOCS.Cmat = 1e-10,
                      epsALS.grandX = 1e-10, epsPOCS.grandX = 1e-10){

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
  sEIG.Cmat <- vector(mode = "list", length = components.Cmat)
  if (sparse.Cmat){
    Cmat.tmp <- Cmat
    masses.Cmat.tmp <- masses.Cmat
    nonzero.relative.ind.list <- vector("list", components.Cmat)
    for (i in 1:components.Cmat) {
      sEIG.Cmat[[i]] <- sparseGEIGEN(
        X = Cmat.tmp, W = masses.Cmat.tmp,
        k = as.integer(components.Cmat),
        init = init.Cmat, seed = seed,
        rds = rds.Cmat[i],
        grp = grp.Cmat,
        orthogonality = orthogonality.Cmat,
        OrthSpace = OrthSpace.Cmat,
        projPriority = projPriority.Cmat,
        itermaxALS = itermaxALS.Cmat, itermaxPOCS = itermaxPOCS.Cmat,
        epsALS = epsALS.Cmat, epsPOCS = epsPOCS.Cmat)

      nonzero.relative.ind.list[[i]] <- which(sEIG.Cmat[[i]]$u[,1] != 0)
      if (length(nonzero.relative.ind.list[[i]]) == 0) stop("Something went wrong !")
      Cmat.tmp <- Cmat.tmp[-nonzero.relative.ind.list[[i]], -nonzero.relative.ind.list[[i]]]
      masses.Cmat.tmp <- masses.Cmat.tmp[-nonzero.relative.ind.list[[i]]]
      if (any(dim(Cmat.tmp) == 0)) stop("Something else went wrong !")
    }
    reference <- 1:nrow(Cmat)
    nonzero.abs.ind.list <- vector("list", length(nonzero.relative.ind.list))
    for (i in 1:components.Cmat) {
      nonzero.abs.ind.list[[i]] <- reference[nonzero.relative.ind.list[[i]]]
      reference <- reference[-nonzero.relative.ind.list[[i]]]
    }
    U <- matrix(0, nrow = nrow(Cmat), ncol = components.Cmat, dimnames = list(rownames(Cmat), paste0("Dim", 1:components.Cmat)))
    U[, 1] <- sEIG.Cmat[[1]]$u[, 1]
    if (components.Cmat > 1) {
      for (i in 2:components.Cmat) {
        U[-nonzero.abs.ind.list[[i - 1]], i] <- sEIG.Cmat[[i]]$u[, 1]
      }
    }
    values <- diag(t(U) %*% Cmat %*% U)
    sEIG.Cmat <- list(values = values,
                      l = values,
                      vectors = U,
                      u = U,
                      rds = rds.Cmat,
                      f = t(t(U * masses.Cmat) * sqrt(values)),
                      iter = NULL,
                      SI = NULL)
  }else{
    ## plain (generalized) eigen
    if (!is.null(masses.Cmat)){
      ## for generalized eigendecomposition (if needed)
      sqrt_W.Cmat <- sqrt(masses.Cmat)
      weighted.Cmat <- t(t(Cmat * sqrt_W.Cmat) * sqrt_W.Cmat)
    }else{
      weighted.Cmat <- Cmat
    }

    eigen.Cmat <- eigen(Cmat, symmetric = TRUE)

    if (!is.null(masses.Cmat)){
      ## for generalized eigendecomposition (if needed)
      eigen.Cmat$u <- eigen.Cmat$vectors
      eigen.Cmat$vectors <- eigen.Cmat$vectors / sqrt_W.Cmat
    }

    sEIG.Cmat <- list(values = eigen.Cmat$values,
                      l = eigen.Cmat$values,
                      vectors = eigen.Cmat$vectors,
                      u = eigen.Cmat$u,
                      rds = rds.Cmat,
                      f = t(t(eigen.Cmat$vectors * masses.Cmat) * sqrt(eigen.Cmat$values)),
                      iter = NULL,
                      SI = NULL)
  }

  ## get alphas

  alpha4grandX <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])
  ### LEFT HERE!
  if (sparse.grandX){
    if (sparseOption == "variable"){
      grpRight = NULL
    }else if (sparseOption == "subtable"){
      grpRight = column.design
    }
  }else{
    rdsLeft.grandX = nrow(data.proc)
    rdsRight.grandX = ncol(data.proc)
  }


  LW <- rep(1/nrow(data.proc), nrow(data.proc))
  RW <- alpha4grandX

  sGSVD.res <- sparseGSVD(X = data.proc, LW = LW, RW = RW, k = components.grandX,
                          init = init.grandX, initLeft = initLeft.grandX, initRight = initRight.grandX, seed = seed,
                          rdsLeft = rdsLeft.grandX, rdsRight = rdsRight.grandX,
                          grpLeft = grpLeft.grandX, grpRight = grpRight.grandX,
                          orthogonality = orthogonality.grandX,
                          OrthSpaceLeft = OrthSpaceLeft.grandX, OrthSpaceRight = OrthSpaceRight.grandX,
                          projPriority = projPriority.grandX,
                          projPriorityLeft = projPriorityLeft.grandX,
                          projPriorityRight = projPriorityLeft.grandX,
                          itermaxALS = itermaxALS.grandX, itermaxPOCS = itermaxPOCS.grandX,
                          epsALS = epsALS.grandX, epsPOCS = epsPOCS.grandX)

  class(sGSVD.res) <- c("sSVD", "sGSVD", "sEIGEN", "sGEIGEN", "sGPCA", "MultiTab", "list")

  res <- spafac.out(sGSVD.res, X = data, LW = LW, RW = RW, tab.idx = tab.idx, column.design = column.design)
  res$X.preproc <- X
  res$alpha <- RW
  res$table.svd <- tab.svd
  res$table.d <- tab.d

  return(res)
}
