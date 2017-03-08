#' compMS2Miner: a package to identify/ visualize unknowns in metabolomic datasets based on MS2 fragmentation data.
#'
#' @description Matches MS1 features to MS2 spectra (.mzXML) files based on a 
#'mass-to-charge and retention time tolerance. Composite spectra and other data
#'can subsequently be visualized during any stage of the compMS2Miner
#'processing workflow. Composite spectra can be denoised, ion signals grouped 
#'and summed, substructure groups identified, common Phase II metabolites
#'predicted and features matched to data bases monoisotopic mass data 
#'and insilico MS2 fragmentation data.
#'The resulting data can then be readily curated by sending to a local or online
#'couchDB database.
#' 
#' @details An example workflow is available in the following vignette:
#' compMS2MinerWorkFlow (source, pdf)
#' 
#' @author WMB Edmands \url{edmandsw@@berkeley.edu}
#' @docType package
#' @name compMS2Miner
NULL
