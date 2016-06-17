#' Metabolite identification methods
#' 
#' @description methods to facilitate metabolite identification including database
#' monoisotopic mass matching, probable annotation filtration, mammalian Phase II
#' metabolite prediction, lipophilicity (logD) prediction, insilico metabolite 
#' fragmentation and metabolite chemical similarity scoring.
#'  
#' @param object. a compMS2 class object
#' @param method. method to use for metabolite identification. See details. 
#' "predSMILES" based on 
#' @param ... option arguments to be passed along.
#' @details Available methods:
#' 
#' 1. monoisotopic mass annotation to data base resources (\code{\link{metID.dbAnnotate}}),
#'    currently available databases include HMDB, DrugBank, T3DB and ReSpect.
#'    possible metabolites electrospray adducts and substructure mass shifts 
#'    are taken into account.
#' 
#' 2. identifies most probable database annotations (\code{\link{metID.dbProb}}), 
#'    taking into account substructure annotations identified by \code{\link{subStructure.Annotate}}.  
#' 
#' 3. Phase II metabolite identification from canonical SMILES currently 
#'    available phase II metabolite prediction types include: acyl-, hydroxl- 
#'    and amine- sulfates and glucuronides and glycine conjugates 
#'    (\code{\link{metID.PredSMILES}}). 
#' 
#' 4. Lipophilicity (LogD) calculation using the ChemAxon command line interface,
#'    this requires ChemAxon software to be installed (In development) (\code{\link{metID.LogDchemAxon}}).
#' 
#' 5. Combinatorial insilico fragment prediction using the command line version
#'    of MetFrag (\code{\link{metID.metFrag}}).
#' 
#' 6. Inter- and intra-feature chemical similarity scoring (In development) (\code{\link{metID.chemSim}}).
#' 
#' 7. Correlation network from a peak table. This function calculates a correlation matrix from the peak areas/ height sample columns and creates a prefuse force directed correlation network that can then be visualized in the \code{\link{compMS2explorer}} application. \code{\link{metID.corrNetwork}}
#' 
#' 8. Spectral similarity network. Inter-spectrum spectral similarity scores (dot product) are calculated. Both fragment ion and precursor - fragment neutral loss
#' pattern similarity scores are calculated and used to identify clusters of spectra with similar fragmentation/neutral loss patterns. A spectral similarity network is then calculated based on a minimum dot product score (minDotProdThresh, default = 0.8).
#' The resulting network can then be visualized in the \code{\link{compMS2explorer}} application. \code{\link{metID.specSimNetwork}}.
#' 
#' @return A compMS2 object with various metabolite identification information.
#' @seealso \code{\link{metID.dbAnnotate}}, \code{\link{metID.dbProb}}, 
#' \code{\link{metID.predSMILES}}, 
#' \code{\link{metID.LogDChemAxon}}, \code{\link{metID.metFrag}}, 
#' \code{\link{metID.chemSim}}, \code{\link{metID.corrNetwork}} 
#' @export
setGeneric("metID", function(object, ...) standardGeneric("metID"))

setMethod("metID", signature = "CompMS2", function(object, method="dbAnnotate", 
                                                    ...) {
  
  method <- match.arg(method, c("dbAnnotate", "dbProb", "predSMILES", 
                                "LogDChemAxon", "metFrag", "CFM", "chemSim", "corrNetwork", 'specSimNetwork'))
  method <- paste("metID", method, sep=".")
  invisible(do.call(method, alist(object, ...)))
}) 