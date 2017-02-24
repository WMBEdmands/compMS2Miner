#' add no MS2 data to compMS2 class object internal to corrNetwork function
#' 
#' @param object a "compMS2" class object.
#' @param specNames character vector of composite spectrum names.
#' @param eicMzRt matrix of EICnos/unique id, mz values and rt values in 3 columns.
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
  if(!is.matrix(eicMzRt)){
    stop('argument eicMzRt must be a matrix') 
  }
network(object) <- list()
emptyList <- vector('list', length(specNames))
names(emptyList) <- specNames
compSpectra(object) <- c(compSpectra(object), emptyList)

metaEleNames <- paste0(rep(specNames, each=3), c('_MS1_EICno', '_MS1_mz', '_MS1_RT'))
eicMzRt <- as.vector(t(eicMzRt))
names(eicMzRt) <- metaEleNames
eicMzRt[grep('_EICno$', metaEleNames)] <- round(eicMzRt[grep('_EICno$', metaEleNames)], 0)
metaDataTmp <- split(eicMzRt, rep(specNames, each=3))
metaData(object) <- c(metaData(object), metaDataTmp)

if(length(DBanno(object)) > 0){
  DBanno(object) <- c(DBanno(object), emptyList)
}
if(length(BestAnno(object)) > 0){
  BestAnno(object) <- c(BestAnno(object), emptyList)
}
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
