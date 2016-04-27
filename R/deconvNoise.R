#' Filter spectral noise from a CompMS2 class object 
#'
#' @param object. a compMS2 class object obtained from the function CompMSset
#' @param method. dynamic noise file "DNF" or fixed maximum intensity "maxInt"
#' @param ... option arguments to be passed along.
#' 
#' @return A compMS2 object with noise filtered composite spectra.
#' @seealso deconvNoise.DNF, deconvNoise.maxInt
#' @export
setGeneric("deconvNoise", function(object, ...) standardGeneric("deconvNoise"))

setMethod("deconvNoise", signature = "CompMS2", function(object, method="DNF", ...) {
  
  method <- match.arg(method, c("DNF","maxInt"))
  method <- paste("deconvNoise", method, sep=".")
  invisible(do.call(method, alist(object, ...)))
})