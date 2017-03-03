#' add no MS2 data to compMS2 class object internal to corrNetwork function
#' 
#' @param object a "compMS2" class object.
#' @param specNames character vector of composite spectrum names.
#' @param eicMzRt data.frame of EICnos/unique id, mz values, rt values and (if applicable)
#' ESI adducts/in-source fragments in 4 columns.
#' @return a "compMS2" class object with noMS2 data added to the appropriate slots.
setGeneric("addNoMS2", function(object, ...) standardGeneric("addNoMS2"))

setMethod("addNoMS2", signature = "compMS2", function(object, specNames=NULL,
                                                      eicMzRt=NULL, ...){
# error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  stopifnot(!is.null(specNames))
  if(!is.character(specNames)){
    stop('argument specNames must be a character vector') 
  }
  stopifnot(!is.null(eicMzRt))
  if(!is.data.frame(eicMzRt)){
    stop('argument eicMzRt must be a data.frame') 
  }
  warning('N.B. You must recalculate any networks after adding features with no MS2 spectra matched.\n', immediate. = TRUE)
network(object) <- list()
emptyList <- vector('list', length(specNames))
names(emptyList) <- specNames
compSpectra(object) <- c(compSpectra(object), emptyList)

metaEleNames <- paste0(rep(specNames, each=4), c('_MS1_EICno', '_MS1_mz', '_MS1_RT', '_MS1_adduct'))
eicMzRt <- as.vector(t(eicMzRt))
names(eicMzRt) <- metaEleNames

metaDataTmp <- split(eicMzRt, rep(specNames, each=4))
# convert to integer and numeric
metaDataTmp <- lapply(metaDataTmp, function(x){
  tmpNames <- names(x)
  names(x) <- NULL
  x <- split(x, 1:4)
  names(x) <- tmpNames
  x[[1]] <- as.integer(x[[1]])
  x[[2]] <- round(as.numeric(x[[2]]), 4)
  x[[3]] <- round(as.numeric(x[[3]]), 3)
  return(x)
})
metaData(object) <- c(metaData(object), metaDataTmp)

if(length(DBanno(object)) > 0){
  DBanno(object) <- c(DBanno(object), emptyList)
}
# if(length(BestAnno(object)) > 0){
#   BestAnno(object) <- c(BestAnno(object), emptyList)
# }
if(nrow(Comments(object)) > 0){
  commentsTmp <- Comments(object)
  emptyDf <- data.frame(matrix('', nrow=length(specNames), ncol=ncol(commentsTmp)), 
                        stringsAsFactors = FALSE)
  colnames(emptyDf) <- colnames(commentsTmp)
  emptyDf$compSpectrum <- specNames
  commentsTmp <- rbind(commentsTmp, emptyDf)
  Comments(object) <- commentsTmp
}
  return(object)
}) # end function
