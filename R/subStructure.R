#' Identify substructures within composite spectra
#' @param object. a compMS2 class object obtained from the function CompMSset
#' @param method. "Annotate" annotation of possible substructure neutral losses/
#' fragments in composite spectra, "prob" identify most probable substructure
#' identification for a composite spectra and "probSummary" summary of probable
#' substructure annotations for each composite spectrum.
#' @param ... option arguments to be passed along.
#' 
#' @return A compMS2 object with substructure annotated composite spectra.
#' @seealso \link{subStructure.Annotate}, \link{subStructure.prob}, \link{subStructure.probSummary}
#' @export
setGeneric("subStructure", function(object, ...) standardGeneric("subStructure"))

setMethod("subStructure", signature = "compMS2", function(object, method="Annotate", 
                                                    ...) {
  
  method <- match.arg(method, c("Annotate","prob","probSummary"))
  method <- paste("subStructure", method, sep=".")
  invisible(do.call(method, alist(object, ...)))
}) 
