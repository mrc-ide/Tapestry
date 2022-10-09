#------------------------------------------------
#' @title Tapestry
#'
#' @description Malaria infections often contain multiple genotypes. When
#'   sequenced together, these produce a complex signal that is a mixture of the
#'   individual genotypes. Like earlier programs like DEploid, Tapestry attempts
#'   to pull these individual genotypes apart by exploiting allele frequency
#'   imbalances within a sample. The difference in Tapestry is that we assume a
#'   model of allele frequencies rather than a reference panel. We also
#'   implement parallel tempered MCMC to ensure good mixing, and produce a
#'   different set of outputs.
#'
#' @docType package
#' @name Tapestry
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib Tapestry, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("Tapestry", libpath)
}
