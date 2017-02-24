#' data.frame of customizable metFrag adduct types and codes
#' 
#' This dataset is the default adduct type table for the \code{\link{metID.metFrag}}
#' function. A custom table can be created following this format as more adduct
#' types are added in future versions of the metFrag command line tool.
#' N.B. adduct names must match those supplied to the \code{\link{adduct2mass}} function internal
#' to the \code{\link{metID.dbAnnotate}} function.
#'
#' \enumerate{
#'  \item adduct adduct name string must match that supplied to adduct2mass.
#'  \item metFragCode the MetFrag command line tool code for the adduct.
#'  \item mode polarity. must be either 'pos' or 'neg'. 
#'  }
#' @docType data
#' @keywords datasets
#' @name metFragAdducts
#' @usage data(metFragAdducts)
#' @source \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/} 
#' @format A data frame with 9 rows and 3 columns
NULL
