#' Sparse Multidimensional Scaling
#'
#' @param DATA
#' @param masses a vector of masses for the rows/columns. Default is set to "NULL" with equal masses of 1/(the number of rows) being used.
#' @param components
#' @param method "distance" or "covariance"
#' @param tol
#' @param init
#' @param seed
#' @param orthogonality
#' @param projPriority
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#' @param rds
#' @param OrthSpace
#'
#' @example
sparseMDS <- function(
  DATA, masses = NULL, components = 0, tol = .Machine$double.eps,
  method="distance",
  init = "svd", seed = NULL,
  grp = NULL,
  rds = rep(1, components),
  orthogonality = "loadings",
  OrthSpace = NULL,
  projPriority = "orth",
  itermaxALS = 1000, itermaxPOCS = 1000,
  epsALS = 1e-10, epsPOCS = 1e-10) {

  ## masses
  if (is.null(masses)){
    M_vec <- rep(1/nrow(DATA), nrow(DATA)) # sum of M_vec = 1
  }else{
    M_vec <- masses
  }

  ## centering matrix
  Cm <- diag(rep(1,dim(DATA)[1])) - (as.matrix(rep(1,dim(DATA)[1])) %*% t(M_vec))

  ## double-centering
  if (method == "covariance"){
    DATA.proc <- 0.5*(Cm %*% DATA %*% t(Cm))
  }else{
    DATA.proc <- -0.5*(Cm %*% DATA %*% t(Cm))
  }

  sGEIG.res <- sparseGEIGEN(X = DATA.proc, W = M_vec, k = components,
                           init = init, seed = seed,
                           rds = rds,
                           grp = grp,
                           orthogonality = orthogonality,
                           OrthSpace = OrthSpace,
                           projPriority = projPriority,
                           itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
                           epsALS = epsALS, epsPOCS = epsPOCS)

  sGEIG.res$input <- list(components = components, init = init,
                         rds = rds, grp = grp,
                         orthogonality = orthogonality, OrthSpace = OrthSpace,
                         projPriority = projPriority,
                         itermaxPOCS = itermaxPOCS, epsPOCS = epsPOCS)

  class(sGEIG.res) <- c("sEIGEN", "sGEIGEN", "sMDS", "list")
  res <- spafac.out(sGEIG.res, X = DATA, LW = M_vec, RW = M_vec)
  res$data <- DATA
  res$data.proc <- DATA.proc
  return(res)

}
