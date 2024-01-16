#' Sparse Multidimensional Scaling
#'
#' Performs sparse multidimensional scaling (MDS) on a given data matrix.
#'
#' @param DATA The data matrix for MDS analysis.
#' @param masses A vector of masses for the rows/columns; defaults to NULL with equal masses of 1/(number of rows).
#' @param components The number of dimensions for the analysis; defaults to 0.
#' @param method Method for analysis, either "distance" or "covariance".
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param init Initialization method; defaults to "svd".
#' @param seed Seed for random number generation; defaults to NULL.
#' @param orthogonality Type of orthogonality constraint; defaults to "loadings".
#' @param projPriority Priority of the projection; defaults to "orth".
#' @param itermaxALS Maximum number of ALS iterations; defaults to 1000.
#' @param itermaxPOCS Maximum number of POCS iterations; defaults to 1000.
#' @param epsALS Convergence criterion for ALS; defaults to 1e-10.
#' @param epsPOCS Convergence criterion for POCS; defaults to 1e-10.
#' @param rds Radii for the analysis; defaults to rep(1, components).
#' @param OrthSpace Orthogonal space for the analysis; defaults to NULL.
#'
#' @return Returns an object containing the results of the sparse MDS analysis, including eigenvalues, eigenvectors, and various metrics relevant to the analysis.
#'
#' @export
#'
#' @examples
#' # Example usage of sparseMDS function
#' # Assuming `data` is a data matrix
#' #' \dontrun{result <- sparseMDS(DATA = data)}
#'
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
