#' Sparse STATIS
#'
#' Performs sparse STATIS analysis on a data matrix with a specific column design.
#'
#' @param X Data matrix with I rows and J columns.
#' @param column.design Design vector indicating the grouping of columns.
#' @param mfa.normalize Logical indicating whether to normalize data for MFA; defaults to TRUE.
#' @param sparse.Cmat Logical indicating whether the C matrix (RV matrix) is sparsified; defaults to FALSE.
#' @param sparse.grandX Logical indicating whether the compromise is sparsified; defaults to TRUE.
#' @param components.Cmat Number of dimensions to extract from the C matrix (RV matrix); defaults to 0.
#' @param components.grandX Number of dimensions to extract from the compromise; defaults to 0.
#' @param sparseOption "variable" or "subtable" indicating the sparsity option.
#' @param center Logical or numeric vector for centering each column of X; passed to \code{\link{scale}}.
#' @param scale Logical or numeric vector for scaling each column of X; passed to \code{\link{scale}}.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param init.Cmat Initialization method for Cmat, and it should be "svd", "eigen", or "rand"; defaults to "eig".
#' @param init.grandX Initialization method for grandX; defaults to "svd".
#' @param initLeft.grandX Initial values for the left side of grandX; defaults to NULL.
#' @param initRight.grandX Initial values for the right side of grandX; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft.Cmat Radii for the left side of Cmat; defaults to rep(1, components.Cmat).
#' @param rdsRight.Cmat Radii for the right side of Cmat; defaults to rep(1, components.Cmat).
#' @param rdsLeft.grandX Radii for the left side of grandX; defaults to rep(1, components.grandX).
#' @param rdsRight.grandX Radii for the right side of grandX; defaults to rep(1, components.grandX).
#' @param grp.Cmat Grouping vector for the Cmat; defaults to NULL.
#' @param grpLeft.grandX Grouping vector for the left side of grandX; defaults to NULL.
#' @param grpRight.grandX Grouping vector for the right side of grandX; defaults to NULL.
#' @param orthogonality.Cmat Type of orthogonality constraint for Cmat; defaults to "loadings".
#' @param OrthSpaceLeft.Cmat Orthogonal space for the left side of Cmat; defaults to NULL.
#' @param OrthSpaceRight.Cmat Orthogonal space for the right side of Cmat; defaults to NULL.
#' @param orthogonality.grandX Type of orthogonality constraint for grandX; defaults to "loadings".
#' @param OrthSpaceLeft.grandX Orthogonal space for the left side of grandX; defaults to NULL.
#' @param OrthSpaceRight.grandX Orthogonal space for the right side of grandX; defaults to NULL.
#' @param projPriority.Cmat Priority of the projection for Cmat; defaults to "orth".
#' @param projPriorityLeft.Cmat Priority of the left projection for Cmat; defaults to projPriority.Cmat.
#' @param projPriorityRight.Cmat Priority of the right projection for Cmat; defaults to projPriority.Cmat.
#' @param projPriority.grandX Priority of the projection for grandX; defaults to "orth".
#' @param projPriorityLeft.grandX Priority of the left projection for grandX; defaults to projPriority.grandX.
#' @param projPriorityRight.Slus Priority of the right projection for grandX; defaults to projPriority.grandX.
#' @param itermaxALS.Cmat Maximum number of ALS iterations for Cmat; defaults to 1000.
#' @param itermaxPOCS.Cmat Maximum number of POCS iterations for Cmat; defaults to 1000.
#' @param itermaxALS.grandX Maximum number of ALS iterations for grandX; defaults to 1000.
#' @param itermaxPOCS.grandX Maximum number of POCS iterations for grandX; defaults to 1000.
#' @param epsALS.Cmat Convergence criterion for ALS for Cmat; defaults to 1e-10.
#' @param epsPOCS.Cmat Convergence criterion for POCS for Cmat; defaults to 1e-10.
#' @param epsALS.grandX Convergence criterion for ALS for grandX; defaults to 1e-10.
#' @param epsPOCS.grandX Convergence criterion for POCS for grandX; defaults to 1e-10.
#'
#' @return Returns an object containing the results of the sparse STATIS analysis, including details of the C matrix (RV matrix), compromise, and various parameters and settings used in the analysis.
#'
#' @export
#'
#' @source Some lines of the function are inspired by the MExPosition package by Derek Beaton and Cherise Chin Fatt.
#' @examples
#' # Example usage of sparseSTATIS function
#' # Assuming `X` is a data matrix and `column.design` is the design vector
#' \dontrun{result <- sparseSTATIS(X = X, column.design = column.design)}
#'
sparseSTATIS <- function(X, column.design, mfa.normalize = TRUE, Cmat.is.RV = TRUE,
                         masses.Cmat = NULL, masses.grandX = NULL,
                      sparse.Cmat = FALSE, sparse.grandX = TRUE,
                      components.Cmat = 2, components.grandX = 2,
                      sparseOption = "variable",
                      center = TRUE, scale = TRUE,
                      tol = .Machine$double.eps,
                      init.Cmat = "svd", init.grandX = "svd",
                      initLeft.grandX = NULL, initRight.grandX = NULL, seed = NULL,
                      rds.Cmat = rep(1, components.Cmat),
                      rdsLeft.grandX = rep(1, components.grandX), rdsRight.grandX = rep(1, components.grandX),
                      grp.Cmat = NULL, grpLeft.grandX = NULL, grpRight.grandX = NULL,
                      orthogonality.Cmat = "loadings",
                      orthogonality.grandX = "loadings",
                      OrthSpace.Cmat = NULL,
                      OrthSpaceLeft.grandX = NULL, OrthSpaceRight.grandX = NULL,
                      projPriority.Cmat = "orth",
                      projPriority.grandX = "orth",
                      projPriorityLeft.grandX = projPriority.grandX,
                      projPriorityRight.grandX = projPriority.grandX,
                      itermaxALS.Cmat = 1000, itermaxPOCS.Cmat = 1000,
                      itermaxALS.grandX = 1000, itermaxPOCS.grandX = 1000,
                      epsALS.Cmat = 1e-10, epsPOCS.Cmat = 1e-10,
                      epsALS.grandX = 1e-10, epsPOCS.grandX = 1e-10){

  if ( !is.matrix(X) ){
    X <- as.matrix(X,rownames.force = TRUE)
  }

  if (sparse.Cmat == TRUE & sparse.grandX == TRUE){
    stop("This function does not do double sparsification.")
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
      data.proc[,tab.idx[1,i]:tab.idx[2,i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]/(tab.d[i]^2) # in a grand table
      data.proclist[[i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]/(tab.d[i]^2) # in a list
    }
  }else{
    alpha.vec <- rep(1, length(column.design))
    data.proc <- data
    for (i in 1:length(count.col)){
      data.proclist[[i]] <- data[,tab.idx[1,i]:tab.idx[2,i]]
    }
  }

  # compute RV matrix an
  Cmat <- GetCmat(data.proclist, RV = Cmat.is.RV, isCube = FALSE)
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
      if (length(nonzero.relative.ind.list[[i]]) == 0) stop("Something went wrong! No element is larger than 0 in the eigenvector of RV.")
      Cmat.tmp <- Cmat.tmp[-nonzero.relative.ind.list[[i]], -nonzero.relative.ind.list[[i]]]
      masses.Cmat.tmp <- masses.Cmat.tmp[-nonzero.relative.ind.list[[i]]]
      if (any(dim(Cmat.tmp) == 0)) stop("Something else went wrong for masses of C because no element is selected according to the sparse RV.")
    }
    # the non-negative eigenvectors
    reference <- 1:nrow(Cmat)
    nonzero.abs.ind.list <- vector("list", length(nonzero.relative.ind.list))
    for (i in 1:components.Cmat) {
      nonzero.abs.ind.list[[i]] <- reference[nonzero.relative.ind.list[[i]]]
      reference <- reference[-nonzero.relative.ind.list[[i]]]
    }
    # create the U matrix to save the sparsified eigenvector for each subspace
    U <- matrix(0, nrow = nrow(Cmat), ncol = components.Cmat, dimnames = list(rownames(Cmat), paste0("Dim", 1:components.Cmat)))
    # for only one subspace
    U[, 1] <- abs(sEIG.Cmat[[1]]$u[, 1]) # flip to positive if they came out negative
    # for more than 1 subspace
    if (components.Cmat > 1) {
      for (i in 2:components.Cmat) {
        U[-nonzero.abs.ind.list[[i - 1]], i] <- abs(sEIG.Cmat[[i]]$u[, 1])
      }
    }
    # create list to save iter, SI, fit ratio, and zero ratio for each subspace
    iter.save <- SI.save <- fitRatio.save <- zeroRatio.save <- list()
    # save output
    for (i in 1:components.Cmat) {
      iter.save[[paste0("Subspace.",i)]] <- sEIG.Cmat[[i]]$iter[1,]
      SI.save[[paste0("Subspace.",i)]] <- sEIG.Cmat[[i]]$SI$SI[1]
      fitRatio.save[[paste0("Subspace.",i)]] <- sEIG.Cmat[[i]]$SI$fitRatio[1]
      zeroRatio.save[[paste0("Subspace.",i)]] <- sEIG.Cmat[[i]]$SI$zeroRatio[1]
    }

    values <- diag(t(U) %*% Cmat %*% U)

    sEIG.Cmat <- list(values = values,
                      l = values,
                      vectors = U,
                      u = U,
                      rds = rds.Cmat,
                      f = t(t(U * masses.Cmat) * sqrt(values)),
                      iter = iter.save,
                      SI = SI.save,
                      fitRatio = fitRatio.save,
                      zeroRatio = zeroRatio.save)
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
                      vectors = eigen.Cmat$vectors,  # (generalized) eigen vectors
                      u = eigen.Cmat$u, # eigen vectors
                      rds = rds.Cmat,
                      f = t(t(eigen.Cmat$vectors * masses.Cmat) * sqrt(eigen.Cmat$values)),
                      iter = NULL,
                      SI = NULL,
                      fitRatio = NULL,
                      zeroRatio = NULL)
  }

  ## get alphas
  if (sparse.Cmat){ # if we sparsified the RV space
    # alphs for difference subspace are stored in different columns
    alpha4grandX <- sweep(sEIG.Cmat$vectors,2,colSums(sEIG.Cmat$vectors),`/`)
  }else{ # if we run a plain eigen on the RV space
    alpha4grandX <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])
  }

  ### Check parameters for sparsifying the grand table
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

  if (sparse.grandX == TRUE){
    if (sparse.Cmat == TRUE){ # if we sparsify both RV and compromise

    }else{ # if we only sparsify the compromise
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

    }
  }else{
    if (sparse.Cmat){ # if we only sparsify the RV (therefore still with the subspaces)

    }else{ # if there is no sparsification

    }
   }


  class(sGSVD.res) <- c("sSVD", "sGSVD", "sEIGEN", "sGEIGEN", "sGPCA", "MultiTab", "list")

  res <- spafac.out(sGSVD.res, X = data, LW = LW, RW = RW, tab.idx = tab.idx, column.design = column.design)
  res$X.preproc <- X
  res$alpha <- RW
  res$table.svd <- tab.svd
  res$table.d <- tab.d

  return(res)
}
