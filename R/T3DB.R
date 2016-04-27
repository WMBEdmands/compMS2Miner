#' T3DB: Toxin and the toxin-target database
#'
#' This dataset contains fields from the T3DB database 
#' The variables are as follows:
#'
#' \itemize{
#'  \item WebAddress web address for database (when combined with the 
#'  Unique_DB_ID entry this will provide a direct link to the online database
#'  entry e.g. Web address : "www.hmdb.ca/metabolites/" +
#'  Unique_DB_ID : "HMDB00001" = \url{www.hmdb.ca/metabolites/HMDB00001})
#'  \item Unique_DB_ID. T3DB unique reference number (T3D0001 -- T3D5000)
#'  \item name. T3DB entry name
#'  \item monoisotopic_weight. Monoisotopic weight of T3DB entry 
#'  (6.032099 -- 3423.581083)
#'  \item SMILES. Canonical SMILES code of T3DB entry
#'  \item molecular_formula. Molecular formula of T3DB entry 
#'  \item cas_number. CAS registry number of T3DB entry
#'  }
#' 
#' @docType data
#' @keywords datasets
#' @name T3DB
#' @usage data(T3DB)
#' @source \url{http://www.T3DB.ca/downloads}
#' @format A data frame with 3528 rows and 7 variables
NULL