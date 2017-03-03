#' subset compMS2 class object using a vector of spectra names
#' 
#' @param object a "compMS2" class object.
#' @param specNames character vector of composite spectrum names.
#' @param corrNetworkNodes logical should all the first correlation network 
#' nodes of the composite spectrum names also be returned? 
#' (default = FALSE).
#' @return a "compMS2" class object with the composite spectra and all metID 
#' information removed. 
#' @export
setGeneric("subsetCompMS2", function(object, ...) standardGeneric("subsetCompMS2"))

setMethod("subsetCompMS2", signature = "compMS2", function(object, 
                                                           specNames=NULL, 
                                                           corrNetworkNodes=FALSE, 
                                                           ...){
# error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  stopifnot(!is.null(specNames))
  if(!is.character(specNames)){
    stop('argument specNames must be a character vector') 
  }
  if(length(network(object)) > 0){
  if(!require(igraph)){
    stop('The igraph package must be installed to use this function.\n')
  }
  }
  if(corrNetworkNodes == TRUE){
    if(length(network(object)$corrNetworkGraph) == 0){
      stop('The function metID.corrNetwork must be run if all connected nodes should be included.')
    }
  }
  matIndxTmp <- match(specNames, names(compSpectra(object))) 
  if(any(is.na(matIndxTmp))){
  stop('The following composite spectrum names do not match:\n', paste0(specNames[is.na(matIndxTmp)], collapse = '\n'), '\nPlease Check and try again.')  
  }
# network(object) <- list()

if(!is.null(network(object)$corrNetworkGraph)){
corrNetTmp <- network(object)$corrNetworkGraph
corrNetIdx <- match(specNames, names(igraph::V(corrNetTmp)))
corrNetIdx <- corrNetIdx[!is.na(corrNetIdx)]
if(corrNetworkNodes){
  # id first neighbours and add to corrNetIdx
  neighSel <- sapply(corrNetIdx, function(x) names(igraph::neighbors(corrNetTmp, x)))
  # add first neighbours to specNames
  specNames <- unique(c(specNames, do.call(c, neighSel)))
  # subset corr network
  corrNetIdx <- match(specNames, names(igraph::V(corrNetTmp)))
  corrNetIdx <- corrNetIdx[!is.na(corrNetIdx)]
  # subset
  corrNetTmp <- igraph::induced_subgraph(corrNetTmp, corrNetIdx)
} 
# subset layout
layoutTmp <- network(object)$corrLayout
corrNetIdx <- match(gsub('.+_', '', names(igraph::V(corrNetTmp))), layoutTmp[, 3])
corrNetIdx <- corrNetIdx[!is.na(corrNetIdx)]
layoutTmp <- layoutTmp[corrNetIdx, , drop=FALSE]
# add back to object
network(object)$corrLayout <- layoutTmp
network(object)$corrNetworkGraph <- corrNetTmp
}
  
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
