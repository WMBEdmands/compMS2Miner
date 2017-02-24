#' summarizes most probable substructure type within all composite spectra
#' 
#' @param object a compMS2 class object
#' @param n number of top substructure types to print.
#' @param minSumRelInt numeric (default = 30)miminum summed relative intensity to consider a probable
#' substructure type identification. 
#' 
#' @return a named numeric vector of frequency of most probable substructure types
#' identified. The most highly ranked probable substructure type for each
#' composite spectra is based on the largest summed relative intensity explained
#' by the characteristic substructure neutral losses and fragments.  
#' @export
setGeneric("subStructure.probSummary", function(object, ...) standardGeneric("subStructure.probSummary"))

setMethod("subStructure.probSummary", signature = "compMS2", function(object, n = 10, 
                                                     minSumRelInt=30) {
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
  if(nrow(subStrAnno(object)) > 0){
    cat("Substructure annotation summary : \n")
    tmp.df <- subStrAnno(object)
    indx.tmp <- duplicated(tmp.df$compSpecName) == F
    noSubStrDet.indx <- tmp.df$SubStrType != "no substructure detected"
    tmp.df$SumRelInt[noSubStrDet.indx == FALSE] <- 0
    aboveMinSumRI <- as.numeric(tmp.df$SumRelInt) > minSumRelInt
    cat(length(which(indx.tmp==TRUE)), "composite spectra \n")
    SubStr.table <- tmp.df$SubStrType[indx.tmp & noSubStrDet.indx & aboveMinSumRI]
    SubStr.table <- sort(table(SubStr.table), decreasing=TRUE)
    cat(sum(SubStr.table), "substructures identified above minimum sum relative intensity of", 
        minSumRelInt, "\n\n") 
    print(SubStr.table[1:ifelse(n > length(SubStr.table), length(SubStr.table), n)])
    return(SubStr.table)
  } else {
    cat("The \"subStructure.prob\" function has not yet been run")
  }
  }
})  
