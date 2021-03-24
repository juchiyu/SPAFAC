#' Create a partition by variable
#'
#' @param X a qualitative variable data frame or matrix
#'
#' @return Returns the partition by variable of a qualitative data set.
#' @export
#'
#' @examples
#'data(cheese)
#'partition_variables(cheese)

partition_variables <- function (Y) {
  Y <- as.data.frame(Y)
  n <- ncol(Y)
  G = list()
  k <- 0

  for (i in 1:n) {
    lev <- levels(as.factor(Y[, i]))
    G[[i]] <- (1+k):(length(lev)+k)
    k <- length(lev)+k
  }
  return(G)
}


