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
  method="distance", Cmat.is.RV = TRUE,
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

  ## eigen decomposition of RV.matrix
  if (sparse.Cmat){ ## NEED SUPSPACE!!
    sEIG.Cmat <- sparseGEIGEN(X = Cmat, W = masses.Cmat, k = as.integer(components.Cmat),
                 init = init.Cmat, seed = seed,
                 rds = rds.Cmat,
                 grp = grp.Cmat,
                 orthogonality = orthogonality.Cmat,
                 OrthSpace = OrthSpace.Cmat,
                 projPriority = projPriority.Cmat,
                 itermaxALS = itermaxALS.Cmat, itermaxPOCS = itermaxPOCS.Cmat,
                 epsALS = epsALS.Cmat, epsPOCS = epsPOCS.Cmat)
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

  alpha4Splus <- sEIG.Cmat$vectors[,1]/sum(sEIG.Cmat$vectors[,1])

  ## compute compromise
  Splus <- ComputeSplus(DATA.proc, alpha4Splus)

  ## Eigen decomposition of Splus
  if (sparse.Splus){
    sGEIG.Splus <- sparseGEIGEN(X = Splus, W = masses.Splus, k = components.Splus,
                              init = init.Splus, seed = seed,
                              rds = rds.Splus,
                              grp = grp.Splus,
                              orthogonality = orthogonality.Splus,
                              OrthSpace = OrthSpace.Splus,
                              projPriority = projPriority.Splus,
                              itermaxALS = itermaxALS.Splus, itermaxPOCS = itermaxPOCS.Splus,
                              epsALS = epsALS.Splus, epsPOCS = epsPOCS.Splus)
  }else{
    ## plain (generalized) eigen
    if (!is.null(masses.Splus)){
      ## for generalized eigendecomposition (if needed)
      sqrt_W.Splus <- sqrt(masses.Splus)
      weighted.Splus <- t(t(Splus * sqrt_W.Splus) * sqrt_W.Splus)
    }else{
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
