#' is_GSVD
#'
#' Tests if the \code{x} object is of class type "GSVD"
#' @details The three primary functions in the \code{GSVD} package produce an inherited (hierarchical) class structure where all of them are of type "GSVD". Those functions are \code{\link{geigen}}, \code{\link{gsvd}}, and \code{\link{gplssvd}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class GSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_GSVD <- function(x){
  inherits(x, "GSVD")
}

#' is_sSVD
#'
#' Tests if the \code{x} object is of class type "sSVD"
#' @details The three primary functions in the \code{sGSVD} package produce an inherited (hierarchical) class structure where all of them are of type "sSVD". Those functions are \code{\link{sparseSVD}} and \code{\link{sparsePLSC}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sGSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_sSVD<- function(x){
  inherits(x, "sSVD")
}

#' is_sGSVD
#'
#' Tests if the \code{x} object is of class type "sGSVD"
#' @details The three primary functions in the \code{sGSVD} package produce an inherited (hierarchical) class structure where all of them are of type "sGSVD". Those functions are \code{\link{sparseGSVD}}, \code{\link{sparseCA}}, \code{\link{sparseMCA}}, \code{\link{sparseMCA}}, \code{\link{sparseMFA}}, \code{\link{sparseDiCA}}, and \code{\link{sPLSCA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sGSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_sGSVD <- function(x){
  inherits(x, "sGSVD")
}

#' is_sEIGEN
#'
#' Tests if the \code{x} object is of class type "sEIGEN"
#' @details The functions in the \code{sGSVD} package that are analyzed by sparse eigendecomposition \code{sparseEIGEN}. For example, \code{sparseMDS} and \code{sparseDiSTATIS}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sEIGEN, FALSE otherwise.
#'
#' @export
is_sEIGEN<- function(x){
  inherits(x, "sEIGEN")
}

#' is_sGEIGEN
#'
#' Tests if the \code{x} object is of class type "sGEIGEN"
#' @details The functions in the \code{sGSVD} package that are analyzed by sparse generalized eigendecomposition \code{sparseGEIGEN}. For example, \code{sparseMDS} and \code{sparseDiSTATIS}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sGEIGEN, FALSE otherwise.
#'
#' @export
is_sGEIGEN<- function(x){
  inherits(x, "sGEIGEN")
}

#' is_sGPCA
#'
#' Tests if the \code{x} object is of class type "sGPCA"
#' @details The three primary functions in the \code{sGPCA} package produce an inherited (hierarchical) class structure where all of them are of type "sGPCA". The function is \code{\link{sparseMFA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sPCA, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_sGPCA <- function(x){
  inherits(x, "sGPCA")
}

#' is_spls
#'
#' Tests if the \code{x} object is of class type "sGSVD"
#' @details The three primary functions in the \code{sGSVD} package produce an inherited (hierarchical) class structure where all of them are of type "spls". Those functions are \code{\link{sparsePLSC}} and \code{\link{sparsePLSCA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sGSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_sPLS <- function(x){
  inherits(x, "sPLS")
}

#' is_discriminant
#'
#' Tests if the \code{x} object is of class type "discriminant"
#' @details The functions in the \code{sGSVD} package that perform discriminant analysis. The function is \code{\link{sparseDiCA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class discriminant, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_discriminant <- function(x){
  inherits(x, "discriminant")
}

#' is_sca
#'
#' Tests if the \code{x} object is of class type "sca"
#' @details The functions in the \code{sGSVD} package that perform discriminant simple correspondence analysis. The function is \code{\link{sparseDiSCA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sca, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_sca <- function(x){
  inherits(x, "SCA")
}


#' Check if an Object is of Class "MCA"
#'
#' Tests if the provided object `x` is of class type "MCA".
#'
#' @details This function is useful for validating whether an object belongs to the class "MCA",
#' which is relevant for functions in the `sGSVD` package, particularly for those performing
#' discriminant multiple correspondence analysis, such as `sparseDiMCA`.
#'
#' @param x object to test.
#' @return A logical value: `TRUE` if the object is of class "MCA", `FALSE` otherwise.
#'
#' @seealso \code{\link{inherits}} for checking object inheritance.
#'
#' @export
is_mca <- function(x) {
  inherits(x, "MCA")
}

#' is_multitab
#'
#' Tests if the \code{x} object is of class type "MultiTab"
#' @details The three primary functions in the \code{sGSVD} package produce an inherited (hierarchical) class structure where all of them are of type "spls". Those functions are \code{\link{sparsePLSC}} and \code{\link{sparsePLSCA}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sGSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_multitab <- function(x){
  inherits(x, "MultiTab")
}

#' is_sMDS
#'
#' Tests if the \code{x} object is of class type "sMDS"
#' @details The functions in the \code{sGSVD} package that are analyzed by sparse multidimensional scaling \code{sparseMDS}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class sMDS, FALSE otherwise.
#'
#' @export
is_sMDS<- function(x){
  inherits(x, "sMDS")
}
