#' Sparse Multiple Factor Analysis
#'
#' Performs sparse multiple factor analysis (MFA) on a given data matrix.
#'
#' @param X Data matrix with I rows and J columns.
#' @param column.design Design vector indicating the grouping of columns.
#' @param components The number of dimensions to extract; defaults to 0.
#' @param sparseOption "variable" or "subtable" indicating the sparsity option.
#' @param center Logical or numeric vector for centering each column of X; passed to \code{\link{scale}}.
#' @param scale Logical or numeric vector for scaling each column of X; passed to \code{\link{scale}}.
#' @param tol Tolerance for the convergence criterion; defaults to \code{.Machine$double.eps}.
#' @param init Initialization method; defaults to "svd".
#' @param initLeft Initial values for the left side; defaults to NULL.
#' @param initRight Initial values for the right side; defaults to NULL.
#' @param seed Seed for random number generation; defaults to NULL.
#' @param rdsLeft Radii for the left side; defaults to rep(1, components).
#' @param rdsRight Radii for the right side; defaults to rep(1, components).
#' @param grpLeft Grouping vector for the left side; defaults to NULL.
#' @param grpRight Grouping vector for the right side; based on sparseOption.
#' @param orthogonality Type of orthogonality constraint; defaults to "loadings".
#' @param OrthSpaceLeft Orthogonal space for the left side; defaults to NULL.
#' @param OrthSpaceRight Orthogonal space for the right side; defaults to NULL.
#' @param projPriority Priority of the projection; defaults to "orth".
#' @param itermaxALS Maximum number of ALS iterations; defaults to 1000.
#' @param itermaxPOCS Maximum number of POCS iterations; defaults to 1000.
#' @param epsALS Convergence criterion for ALS; defaults to 1e-10.
#' @param epsPOCS Convergence criterion for POCS; defaults to 1e-10.
#'
#' @return Returns an object containing the results of the sparse MFA analysis. This object includes the results of sparse generalized singular value decomposition (sparseGSVD) and additional data specific to MFA.
#'
#' @export
#'
#' @source Some lines of the function are inspired by the MExPosition package by Derek Beaton and Cherise Chin Fatt.
#' @examples
#' # Example usage of sparseMFA function
#' # Assuming `data` is a data matrix and `column.design` is the design vector
#' \dontrun{result <- sparseMFA(X = data, column.design = column.design)}
#'
sparseMFA <- function(X, column.design, components = 0,
                      sparseOption = "variable",
                       center = TRUE, scale = TRUE, mfa.scale = TRUE,
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

  if ( !is.matrix(X) ){
    X <- as.matrix(X,rownames.force = TRUE)
  }

  if (scale == "ss1"){
    data <- projL2(scale(X, center = center, scale = FALSE))$x  # center and scale to ss = 1 before MFA-normalization
  }else{
    data <- scale(X, center = center, scale = scale) # center and scale before MFA-normalization
  }



  if (mfa.scale){
    alpha.vec <- vector(length = length(column.design))
    count.col <- table(column.design)
    ncol.tab <- c(0, count.col)
    tab.idx <- matrix(nrow = 2, ncol = length(count.col), dimnames = list(c("from", "to"), names(count.col)))
    tab.svd <- list()
    tab.d <- vector(length = length(count.col))
    to_total = 0
    from_total = 1
    for (i in 1:length(count.col)) {
      from = ncol.tab[i] + from_total
      to =  ncol.tab[i+1] + to_total
      to_total = to
      from_total = from
      tab.svd[[i]] <- svd(data[, from:to])
      tab.d[i] <- tab.svd[[i]]$d[1] # normalize the center-and-scaled tables
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
                         epsALS = epsALS, epsPOCS = epsPOCS
)

  class(sGSVD.res) <- c("sSVD", "sGSVD", "sGPCA", "MultiTab", "list")

  res <- spafac.out(sGSVD.res, X = data, LW = LW, RW = RW, tab.idx = tab.idx, column.design = column.design)
  res$X.preproc <- X
  res$alpha <- RW
  res$table.svd <- tab.svd
  res$table.d <- tab.d

  return(res)
}
