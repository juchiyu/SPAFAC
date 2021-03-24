


#' Create a partition from a variable
#'
#' @param vect the factor vector used to create the partition. The vector length must be equal to the number of observation or of variables.
#' @param dta the complete data set used in the case of the variables (columns) partition.
#' @param choice Indicate if you want a partition on the observations ("obs") or on the variables ("var") of the data set (default obs).
#'
#' @return Returns a list with the index of each partition levels.
#' @export
#'
#' @examples
#' data(cheese)
#' G <- partition_select (cheese [,3])
#' vect <- sample(c("a", "b"), size = ncol(cheese), replace = T)
#' G2 <- partition_select (vect, dta = cheese, choice = "var")

partition_select <- function (vect, dta = NULL, choice = "obs") {

  lev <- levels(as.factor(vect))
  G <- lapply(lev, function(i) {G <- as.vector(which(vect==i))})
  names(G) <- lev

  if (choice == "var") {
    G2 <- list()
    Gvar <- partition_variables(dta)
    for (i in 1:length(lev)) {
      G2[[i]] <- do.call(c, sapply(G[[i]], function(j) {Gvar[[j]]}))
    }
    return(G2)
  }else{return(G)}

}
