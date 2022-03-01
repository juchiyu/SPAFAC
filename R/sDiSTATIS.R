#' Sparse DiSTATIS
#'
#' @param DATA a cube of distance matrix
#' @param method "distance" or "covariance"
#' @param masses.Cmat a vector of masses for the rows/columns of the rv matric. Default is set to "NULL" with all 1s.
#' @param masses.Splus a vector of masses for the rows/columns of each table in the cube. Default is set to "NULL" with equal masses of 1/(the number of rows) being used.
#' @param sparse.Cmat
#' @param sparse.Splus
#' @param components.Cmat
#' @param components.Splus
#' @param init.Cmat
#' @param init.Splus
#' @param seed
#' @param grp.Cmat
#' @param grp.Splus
#' @param rds.Cmat
#' @param rds.Splus
#' @param orthogonality.Cmat
#' @param OrthSpace.Cmat
#' @param orthogonality.Splus
#' @param OrthSpace.Splus
#' @param projPriority.Cmat
#' @param projPriority.Splus
#' @param itermaxALS.Cmat
#' @param itermaxPOCS.Cmat
#' @param itermaxALS.Splus
#' @param itermaxPOCS.Splus
#' @param epsALS.Cmat
#' @param epsPOCS.Cmat
#' @param epsALS.Splus
#' @param epsPOCS.Splus
#' @param tol
#'
#' @example
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
        apply(DATA, 3, distatis.preproc, method = method, masses = masses.Splus)
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
      stop("The length of `masses.Cmat` does not equal the number of table (i.e., the third dimension of `DATA`).")
    }else{
      masses.Cmat <- masses.Cmat
    }
  }

  ## eigen decomposition of RV.matrix
  sEIG.Cmat <- vector(mode = "list", length = components.Cmat)
  if (sparse.Cmat) { ## NEED SUPSPACE!!

    Cmat.tmp <- Cmat
    nonzero.relative.ind.list <- vector("list", components.Cmat)
    for (i in 1:components.Cmat) {
      sEIG.Cmat[[i]] <- sparseGEIGEN(
        X = Cmat.tmp, W = masses.Cmat,
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
      if (any(dim(Cmat.tmp) == 0)) stop("Something else went wrong !")
    }
    reference <- 1:nrow(Cmat)
    nonzero.abs.ind.list <- vector("list", length(nonzero.relative.ind.list))
    for (i in 1:components.Cmat) {
      nonzero.abs.ind.list[[i]] <- reference[nonzero.relative.ind.list[[i]]]
      reference <- reference[-nonzero.relative.ind.list[[i]]]
    }
    U <- matrix(0, nrow = nrow(Cmat), ncol = components.Cmat, dimnames = list(rownames(Cmat), paste0("Dim", 1:components.Cmat)))
    U[, 1] <- sEIG.Cmat[[i]]$u[, 1]
    if (components.Cmat > 1) {
      for (i in 2:components.Cmat) {
        U[nonzero.abs.ind.list[[i - 1]], i] <- sEIG.Cmat[[i]]
      }
    }
    values <- drop(t(U) %*% Cmat %*% U)
    sEIG.Cmat <- list(values = values,
                      l = values,
                      vectors = U,
                      u = U,
                      rds = rds.Cmat,
                      f = t(t(U * masses.Cmat) * sqrt(values)),
                      iter = NULL,
                      SI = NULL)
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
  }

  ## get alphas

  alpha4Splus <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])

  ## compute compromise
  Splus <- ComputeSplus(DATA.proc, alpha4Splus)
  sGEIG.Splus <- vector(mode = "list", length = components.Cmat)
  ## Eigen decomposition of Splus
  if (sparse.Splus){
    for (i in 1:components.Cmat) {
      sGEIG.Splus[[i]] <- sparseGEIGEN(X = Splus, W = masses.Splus, k = components.Splus,
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
  res <- spafac.out(sGEIG.res, X = DATA, LW = M_vec, RW = M_vec)
  res$data <- DATA
  res$data.proc <- DATA.proc
  return(res)

}



it1 <- list(v = c(1.2, 6.3), i = c(4, 2), Unknown = c(4, 2)) # 123456
it2 <- list(v = c(0.5, 2.2), i = c(4, 2), Unknown = c(6, 3)) # 1356
it3 <- list(v = c(0.1, 0.2), i = c(1, 2), Unknown = c(1, 5)) # 15
x1 <- x2 <- x3 <- rep(0, 6)
x1[it1$i] <- it1$v
x2[-it1$i][it2$i] <- it2$v
x3[-it1$i][-it2$i][it3$i] <- it3$v
cbind(x1, x2, x3)
# [1] 0.1 6.3 2.2 1.2 0.2 0.5
result <- cbind(
  c(0, 6.3, 0, 1.2, 0, 0),
  c(0, 0, 2.2, 0, 0, 0.5),
  c(0.1, 0, 0, 0, 0.2, 0))

itlist <- list(it1, it2, it3)
reference <- 1:6
unknown <- vector("list", 3)
for (i in 1:3) {
  unknown[[i]] <- reference[itlist[[i]]$i]
  reference <- reference[-itlist[[i]]$i]
}


