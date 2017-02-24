#' subset compMS2 class object using a vector of spectra names
#' 
#' @param object a "compMS2" class object.
#' @param specNames character vector of composite spectrum names.
#' @return a "compMS2" class object with the composite spectra and all metID 
#' information removed. Any correlation or spectral similarity networks will 
#' have to be recalculated.
#' @export
setGeneric("subsetCompMS2", function(object, ...) standardGeneric("subsetCompMS2"))

setMethod("subsetCompMS2", signature = "compMS2", function(object, 
                                                           specNames=NULL, ...){
# error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  stopifnot(!is.null(specNames))
  if(!is.character(specNames)){
    stop('argument specNames must be a character vector') 
  }
  matIndxTmp <- specNames %in% names(compSpectra(object)) 
  if(any(matIndxTmp == FALSE)){
  stop('The following composite spectrum names do not match:\n', paste0(specNames[matIndxTmp == FALSE], collapse = '\n'), '\nPlease Check and try again.')  
  }
network(object) <- list()
compSpectra(object) <- compSpectra(object)[specNames]
metaData(object) <- metaData(object)[specNames]
if(length(DBanno(object)) > 0){
  DBanno(object) <- DBanno(object)[specNames]
}
if(length(BestAnno(object)) > 0){
  BestAnno(object) <- BestAnno(object)[specNames]
}
if(nrow(subStrAnno(object)) > 0){
  subStrAnno(object) <- subStrAnno(object)[subStrAnno(object)$compSpecName %in% specNames, , drop=FALSE]
}
if(length(object@spectralDB) > 0){
  object@spectralDB <- object@spectralDB[specNames]
}
if(length(object@inSilico) > 0){
  if(!is.null(object@inSilico$MetFrag)){
    object@inSilico$MetFrag <- object@inSilico$MetFrag[specNames]  
  }
  if(!is.null(object@inSilico$CFM)){
    object@inSilico$CFM <- object@inSilico$CFM[specNames]  
  }
}
if(nrow(Comments(object)) > 0){
  Comments(object) <- Comments(object)[Comments(object)$compSpectrum %in% specNames, , drop=FALSE]
}
  return(object)
}) # end function
