#' Sparse DiSTATIS
#'
#' Performs Sparse DiSTATIS analysis on a cube of distance matrices.
#'
#' @param DATA A cube of distance matrices for analysis.
#' @param method Method for analysis, either "distance" or "covariance".
#' @param masses.Cmat A vector of masses for the rows/columns of the RV matrix; defaults to NULL with all 1s.
#' @param masses.Splus A vector of masses for the rows/columns of each table in the cube; defaults to NULL with equal masses of 1/(number of rows).
#' @param sparse.Cmat Logical indicating whether to apply sparsity on RV matrix; defaults to FALSE.
#' @param sparse.Splus Logical indicating whether to apply sparsity on each table; defaults to TRUE.
#' @param components.Cmat Number of components for RV matrix; defaults to 0.
#' @param components.Splus Number of components for each table; defaults to 0.
#' @param init.Cmat Initialization method for RV matrix; defaults to "svd".
#' @param init.Splus Initialization method for each table; defaults to "svd".
#' @param seed Seed for random number generation; defaults to NULL.
#' @param grp.Cmat Grouping vector for the RV matrix; defaults to NULL.
#' @param grp.Splus Grouping vector for each table; defaults to NULL.
#' @param rds.Cmat Radii for the RV matrix; defaults to rep(1, components.Cmat).
#' @param rds.Splus Radii for each table; defaults to rep(1, components.Splus).
#' @param orthogonality.Cmat Type of orthogonality constraint for RV matrix; defaults to "loadings".
#' @param OrthSpace.Cmat Orthogonal space for the RV matrix; defaults to NULL.
#' @param orthogonality.Splus Type of orthogonality constraint for each table; defaults to "loadings".
#' @param OrthSpace.Splus Orthogonal space for each table; defaults to NULL.
#' @param projPriority.Cmat Priority of the projection for RV matrix; defaults to "orth".
#' @param projPriority.Splus Priority of the projection for each table; defaults to "orth".
#' @param itermaxALS.Cmat Maximum number of ALS iterations for RV matrix; defaults to 1000.
#' @param itermaxPOCS.Cmat Maximum number of POCS iterations for RV matrix; defaults to 1000.
#' @param itermaxALS.Splus Maximum number of ALS iterations for each table; defaults to 1000.
#' @param itermaxPOCS.Splus Maximum number of POCS iterations for each table; defaults to 1000.
#' @param epsALS.Cmat Convergence criterion for ALS for RV matrix; defaults to 1e-10.
#' @param epsPOCS.Cmat Convergence criterion for POCS for RV matrix; defaults to 1e-10.
#' @param epsALS.Splus Convergence criterion for ALS for each table; defaults to 1e-10.
#' @param epsPOCS.Splus Convergence criterion for POCS for each table; defaults to 1e-10.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#'
#' @return Returns an object containing the results of the Sparse DiSTATIS analysis, including results for both the RV matrix and each individual table, along with various parameters and settings used in the analysis.
#'
#' @export
#'
#' @examples
#' # Example usage of sparseSTATIS function
#' # Assuming `X` is a data matrix and `column.design` is the design vector
#' \dontrun{result <- sparseDiSTATIS(X = X, column.design = column.design)}
#'
sparseDiSTATIS <- function(
  DATA,
  method = "distance", Cmat.is.RV = TRUE,
  masses.Cmat = NULL, masses.Splus = NULL,
  sparse.Cmat = FALSE, sparse.Splus = TRUE,
  components.Cmat = 0, components.Splus = 0,
  init.Cmat = "svd", init.Splus = "svd",
  seed = NULL,
  grp.Cmat = NULL, grp.Splus = NULL,
  rds.Cmat = rep(1, components.Cmat), rds.Splus = rep(1, components.Splus),
  orthogonality.Cmat = "loadings", OrthSpace.Cmat = NULL,
  orthogonality.Splus = "loadings", OrthSpace.Splus = NULL,
  projPriority.Cmat = "orth", projPriority.Splus = "orth",
  itermaxALS.Cmat = 1000, itermaxPOCS.Cmat = 1000,
  itermaxALS.Splus = 1000, itermaxPOCS.Splus = 1000,
  epsALS.Cmat = 1e-10, epsPOCS.Cmat = 1e-10,
  epsALS.Splus = 1e-10, epsPOCS.Splus = 1e-10,
  tol = .Machine$double.eps) {

  ## preprocessing
  DATA.proc <- array(
    as.numeric(
      unlist(
        apply(DATA, 3, distatis.preproc, method = method, masses = masses.Splus, simplify = FALSE)
      )
    ), dim = dim(DATA), dimnames = dimnames(DATA))

  ## compute RV matrix and get it's masses
  Cmat <- GetCmat(DATA.proc, RV = Cmat.is.RV)
  ### masses for the RV matrix
  if (is.null(masses.Cmat)) {
    masses.Cmat <- rep(1, nrow(Cmat))
  } else {
    ## check length
    if (length(masses.Cmat) != nrow(Cmat)) {
      stop("The length of `masses.Cmat` does not equal the number of tables (i.e., the third dimension of `DATA`).")
    } else {
      masses.Cmat <- masses.Cmat
    }
  }

  ### masses for the Splus matrix
  if (is.null(masses.Splus)) {
    masses.Splus <- rep(1, nrow(DATA))
  } else {
    ## check length
    if (length(masses.Splus) != nrow(DATA)) {
      stop("The length of `masses.Splus` does not equal the number of subjects (i.e., the frist dimension of `DATA`).")
    }else{
      masses.Splus <- masses.Splus
    }
  }
  ## eigen decomposition of RV.matrix
  sEIG.Cmat <- vector(mode = "list", length = components.Cmat)
  if (sparse.Cmat) {

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
    ## get alphas
    alpha4Splus <- lapply(as.list(data.frame(sEIG.Cmat$vectors[, 1:components.Cmat])), function(x) x/sum(x))

  } else {
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
    ## get alphas
    alpha4Splus <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])
  }

  ### move sparseGEIGEN into the if-else

  ## compute compromise
  Splus <- ComputeSplus(DATA.proc, alpha4Splus)
  sGEIG.Splus <- vector(mode = "list", length = components.Cmat)
  ## Eigen decomposition of Splus
  if (sparse.Splus){

    ## compute compromise
    Splus <- lapply(alpha4Splus, function(alpha) ComputeSplus(DATA.proc, alpha))
    sGEIG.Splus <- vector(mode = "list", length = components.Cmat)

    for (i in 1:components.Cmat) {
      sGEIG.Splus[[i]] <- sparseGEIGEN(X = Splus[[i]],
                                       W = masses.Splus, k = components.Splus,
                                       init = init.Splus, seed = seed,
                                       rds = rds.Splus,
                                       grp = grp.Splus,
                                       orthogonality = orthogonality.Splus,
                                       OrthSpace = OrthSpace.Splus,
                                       projPriority = projPriority.Splus,
                                       itermaxALS = itermaxALS.Splus, itermaxPOCS = itermaxPOCS.Splus,
                                       epsALS = epsALS.Splus, epsPOCS = epsPOCS.Splus)
    }
  } else {

    ## BELOW, we need a LOOP
    alpha4Splus <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])

    ## compute compromise
    Splus <- ComputeSplus(DATA.proc, alpha4Splus)
    sGEIG.Splus <- vector(mode = "list", length = components.Cmat)

    ## plain (generalized) eigen
    if (!is.null(masses.Splus)){
      ## for generalized eigendecomposition (if needed)
      sqrt_W.Splus <- sqrt(masses.Splus)
      weighted.Splus <- t(t(Splus * sqrt_W.Splus) * sqrt_W.Splus)
    } else {
      weighted.Splus <- Splus
    }

    eigen.Splus <- eigen(Splus, symmetric = TRUE)

    if (!is.null(masses.Splus)){
      ## for generalized eigendecomposition (if needed)
      eigen.Splus$u <- eigen.Splus$vectors
      eigen.Splus$vectors <- eigen.Splus$vectors / sqrt_W.Splus
    }

    sGEIG.Splus <- list(values = eigen.Splus$values,
                      l = eigen.Splus$values,
                      vectors = eigen.Splus$vectors,
                      u = eigen.Splus$u,
                      rds = rds.Splus,
                      f = t(t(eigen.Splus$vectors * masses.Splus) * sqrt(eigen.Splus$values)),
                      iter = NULL,
                      SI = NULL)
  }


  sGEIG.res <- list(
    res4Cmat <- sEIG.Cmat,
    res4Cmat$Cmat <- Cmat,
    res4Splus <- sGEIG.Splus,
    res4Splus$Splus <- Splus,
    res4Splus$alpha <- alpha4Splus,
    param <- list(components.Cmat = components.Cmat, init.Cmat = init.Cmat,
                  rds.Cmat = rds.Cmat, grp.Cmat = grp.Cmat,
                  orthogonality.Cmat = orthogonality.Cmat, OrthSpace.Cmat = OrthSpace.Cmat,
                  projPriority.Cmat = projPriority.Cmat,
                  itermaxPOCS.Cmat = itermaxPOCS.Cmat, epsPOCS.Cmat = epsPOCS.Cmat,
                  components.Splus = components.Splus, init.Splus = init.Splus,
                  rds.Splus = rds.Splus, grp.Splus = grp.Splus,
                  orthogonality.Splus = orthogonality.Splus, OrthSpace.Splus = OrthSpace.Splus,
                  projPriority.Splus = projPriority.Splus,
                  itermaxPOCS.Splus = itermaxPOCS.Splus, epsPOCS.Splus = epsPOCS.Splus)
  )

  class(sGEIG.res) <- c("sEIGEN", "sGEIGEN", "MultiTab", "list")
  res <- spafac.out(sGEIG.res, X = DATA, masses.Cmat =  masses.Cmat,
                    LW = masses.Splus, LM = masses.Splus)
  res$data <- DATA
  res$data.proc <- DATA.proc
  return(res)

}

# it1 <- list(v = c(1.2, 6.3), i = c(4, 2), Unknown = c(4, 2)) # 123456
# it2 <- list(v = c(0.5, 2.2), i = c(4, 2), Unknown = c(6, 3)) # 1356
# it3 <- list(v = c(0.1, 0.2), i = c(1, 2), Unknown = c(1, 5)) # 15
# x1 <- x2 <- x3 <- rep(0, 6)
# x1[it1$i] <- it1$v
# x2[-it1$i][it2$i] <- it2$v
# x3[-it1$i][-it2$i][it3$i] <- it3$v
# cbind(x1, x2, x3)
# # [1] 0.1 6.3 2.2 1.2 0.2 0.5
# result <- cbind(
#   c(0, 6.3, 0, 1.2, 0, 0),
#   c(0, 0, 2.2, 0, 0, 0.5),
#   c(0.1, 0, 0, 0, 0.2, 0))
#
# itlist <- list(it1, it2, it3)
# reference <- 1:6
# unknown <- vector("list", 3)
# for (i in 1:3) {
#   unknown[[i]] <- reference[itlist[[i]]$i]
#   reference <- reference[-itlist[[i]]$i]
# }
#
# # Clean workspace
# rm(list = ls()) ; graphics.off()
#
# # Load libraries
# library(devtools)
# library(DistatisR)
#
# load_all(".")
# load_all("../sGSVD/")
#
# # Load data
# data(DistAlgo)
# DistAlgo[3,1,1] <- DistAlgo[1,3,1]
# # Regular DiStatis
# res.reg <- distatis(DistAlgo)
#
# # Sparse DiStatis
# res.sp <- sparseDiSTATIS(DistAlgo, sparse.Cmat = TRUE, sparse.Splus = TRUE, components.Cmat = 2L, components.Splus = 2L)
#
#
# A <- res.reg$res4Splus$Splus
# A - t(A)

