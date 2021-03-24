#' Disjonctif table
#'
#' @param X A data frame or matrix with factors
#'
#' @return Returns the disjonctif table
#' @export
#'
#' @examples
#' library(FactoMineR)
#' data(tea)
#' tea <- tea [, 1:8]
#'
#' tab_disjonctif(tea)
#'
tab_disjonctif <- function (X) {
  X <- as.matrix(X)
  m <- nrow(X)
  n <- ncol(X)

  X <- as.data.frame(X)
  newdata <- lapply (1:n, function (i) {
    lev <- levels(as.factor(X[, i]))

    stmp <- sapply (1:m, function (j) {
      as.numeric(X[j, i] == lev) #the level equal to the row will get 1 and the others 0
    })
    if (length(lev)==1) {
      s <- as.matrix(stmp)
    }else{
      s <- t(stmp)
    }
    colnames(s) <- paste(colnames(X)[i], lev, sep = "_")
    s
  })

  newdata <- do.call(cbind, newdata)
  return(newdata)
}
