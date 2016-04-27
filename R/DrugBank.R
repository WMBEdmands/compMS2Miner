#' DrugBank Database 
#'
#' This dataset contains fields from the DrugBank database 
#' The variables are as follows:
#'
#' \itemize{
#'  \item WebAddress web address for database (when combined with the 
#'  Unique_DB_ID entry this will provide a direct link to the online database
#'  entry e.g. Web address : "www.HMDB.ca/metabolites/" +
#'  Unique_DB_ID : "HMDB00001" = \url{www.HMDB.ca/metabolites/HMDB00001})
#'  \item Unique_DB_ID. DrugBank unique reference number (DB00091 -- DB09028)
#'  \item name. DrugBank entry name
#'  \item monoisotopic_weight. Monoisotopic weight of DrugBank entry 
#'  (7.016004 -- 6176.017855)
#'  \item SMILES. Canonical SMILES code of DrugBank entry
#'  \item molecular_formula. Molecular formula of DrugBank entry 
#'  \item group. drug type group see unique(DrugBank$group)
#'  \item cas_number. CAS registry number of DrugBank entry
#'  }
#' 
#' @docType data
#' @keywords datasets
#' @name DrugBank
#' @usage data(DrugBank)
#' @source \url{http://www.drugbank.ca/downloads}
#' @format A data frame with 6573 rows and 8 variables
NULL