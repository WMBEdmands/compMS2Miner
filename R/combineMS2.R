#' combineMS2 ions either within and/ or between composite spectra 
#' 
#' @param object. a compMS2 class object obtained from the function CompMSset
#' @param method. "Ions" intra-spectrum ions grouping or "Spectra" inter-spectra
#' ion grouping, composite spectra from multiple MS2 files matched to the same
#' MS1 feature are matched.
#' @param ... option arguments to be passed along.
#' 
#' @return A compMS2 object with intra-composite spectrum or inter-composite 
#' spectra, grouped, signal summed and group mass-to-charge averaged spectra. 
#' @seealso \code{\link{combineMS2.Ions}}, \code{\link{combineMS2.Spectra}}
#' @export
setGeneric("combineMS2", function(object, ...) standardGeneric("combineMS2"))

setMethod("combineMS2", signature = "CompMS2", function(object, method="DNF", ...) {
  
  method <- match.arg(method, c("Ions","Spectra"))
  method <- paste("combineMS2", method, sep=".")
  invisible(do.call(method, alist(object, ...)))
})