#' Interfeature correlation
#' @param sampleIDstring Unique text string to identify sample columns for correlation matrix calculation 
setGeneric("intFeatCorr", function(object, ...) standardGeneric("intFeatCorr"))

setMethod("intFeatCorr", signature = "CompMS2", function(object,
                                                         sampleIDstring = NULL){
  
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (is.null(sampleIDstring)){
    stop("argument sampleIDstring is missing with no default")
  } else {
  }
})
