#' Identifies probable substructure type 
#' 
#' @description Identifies probable substructure type based on the summed
#' relative intensites and therefore proportion of total composite spectrum 
#' intensity explained.
#' @param object compMS2 object 
#' @param minSumRelInt numeric (default = 70) miminum summed relative intensity to consider a probable
#' substructure type identification. If above this minimum summed relative intensity
#' then the most probable substructure type will be added to the compound_class 
#' column of the Comments table in the compMS2 object with a note stating it was
#' identified using this function e.g. sulfate (subStructure.prob). This provides
#' a means to automatically annotate the Comments table in a first-pass 
#' metabolite identification workflow. 
#' @return a data.frame of probable substructure annotations for each composite
#' spectrum, ranked by the sum of the relative intensities for that substructure
#' type 
#' @export
setGeneric("subStructure.prob", function(object, ...) standardGeneric("subStructure.prob"))

setMethod("subStructure.prob", signature = "compMS2", function(object, 
                                                               minSumRelInt=70){
  
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if ("Frag.ID.type" %in% colnames(compSpectra(object)[[1]])){
    
    message("Identifying likely substructure type...")
    flush.console()
    
    names_CompSpec <- names(compSpectra(object))
    substrBestIds <- lapply(c(1:length(compSpectra(object))), function(x){
      spec.df <- compSpectra(object)[[x]]
      if(nrow(spec.df) > 1){
      spec.df <- as.data.frame(apply(spec.df, 2, function(y) gsub('noID|noID;|;noID', '', y)), stringsAsFactors = FALSE)
      } else {
        spec.df <- t(as.data.frame(apply(spec.df, 2, function(y) gsub('noID|noID;|;noID', '', y)), stringsAsFactors = FALSE))
        row.names(spec.df) <- NULL
      }
      arrIndxAnno <- which(spec.df[, c("Frag.ID.type", "interfrag.loss.type", 
                                        "Neutral.loss.type")] != '', arr.ind = TRUE)
      subStrTypes <- data.frame(compSpecName = names_CompSpec[x],
                                SubStrType = "no substructure detected", 
                                SumRelInt = "no substructure detected", 
                                Freq = "no substructure detected",
                                nPeaks = "no substructure detected",
                                stringsAsFactors = FALSE)
      if(is.matrix(arrIndxAnno)){
      if(nrow(arrIndxAnno) > 0){
      subStrTypes <- spec.df[, c("Frag.ID.type", "interfrag.loss.type", 
                                 "Neutral.loss.type")]
      subStrTypes <- subStrTypes[arrIndxAnno]
      names(subStrTypes) <- spec.df$Rel_Intensity[arrIndxAnno[, 1]]
      relIntOrder <- order(spec.df$Rel_Intensity[arrIndxAnno[, 1]], decreasing = TRUE)
      fragNLtypes <- spec.df[, c("Frag.ID", "interfrag.loss", 
                                 "Neutral.loss")]
      fragNLnames <- colnames(fragNLtypes)[arrIndxAnno[, 2]]
      fragNLtypes <- fragNLtypes[arrIndxAnno]
      fragNLtypes <- paste0(fragNLnames, '_', fragNLtypes)
      # order rel int
      subStrTypes <- subStrTypes[relIntOrder]
      fragNLtypes <- fragNLtypes[relIntOrder]
      duplIndx <- duplicated(fragNLtypes) == F
      subStrTypes <- subStrTypes[duplIndx]
      # remove duplicates
      # subStrTypes <- paste0(subStrTypes, ';;', 
      #                              gsub("[A-Za-z]|\\.", "", names(subStrTypes)))
      # duplIndx <- duplicated(subStrTypes) == F
      # subStrTypes <- subStrTypes[duplIndx]
      # fragNLTypes <- fragNLTypes[duplIndx]
      # # vector intensities
      # names(subStrTypes) <- spec.df$Rel_Intensity[as.numeric(gsub(".+;;", 
      #                                           "", subStrTypes))]
      # unlist any multiple IDs
      subStrTypes <- unlist(strsplit(subStrTypes, ";"))
      # remove last character string split
      names(subStrTypes) <- substring(names(subStrTypes), 1, 
                                      ifelse(as.numeric(names(subStrTypes)) > 1000,
                                             nchar(names(subStrTypes))- 1, nchar(names(subStrTypes))))
#        <- Ints
      # sum rel intensities
      subStrs.tmp <- data.frame(table(subStrTypes), stringsAsFactors = FALSE)
      subStrTypes <- tapply(as.numeric(names(subStrTypes)), subStrTypes, sum)
           # subStrTypes <- subStrTypes / sum(spec.df$intensity) 
      # summary results
      subStrTypes <- data.frame(compSpecName = rep(names_CompSpec[x], 
                                                   length(subStrTypes)),
                                SubStrType = names(subStrTypes), 
                                SumRelInt = subStrTypes,
                                stringsAsFactors = FALSE)
      subStrTypes <- merge(subStrTypes, subStrs.tmp, by.x = "SubStrType", 
                           by.y = "subStrTypes")
      subStrTypes <- subStrTypes[order(subStrTypes$SumRelInt,
                                       decreasing = TRUE), , drop=FALSE]
      subStrTypes$nPeaks <- length(spec.df$intensity)
      } 
    }
      return(subStrTypes)
    })
    
  } else {
    stop("subStructure.Annotate function has not yet been applied")
  }
     
    subStrAnno(object) <- do.call("rbind", substrBestIds)
    # add any compound classes to the Comments table
    subStrTypes <- subStrAnno(object)
    # subset
    subStrTypes <- subStrTypes[subStrTypes$SumRelInt != "no substructure detected", , drop=FALSE]
    subStrTypes <- subStrTypes[as.numeric(subStrTypes$SumRelInt) >= minSumRelInt, , drop=FALSE]
    subStrTypes <- subStrTypes[duplicated(subStrTypes$compSpecName) == FALSE, , drop=FALSE]
    metIDcomments <- Comments(object)
    indxTmp <- match(metIDcomments$compSpectrum, subStrTypes$compSpecName)
    metIDcomments$compound_class[!is.na(indxTmp)] <- paste0(subStrTypes$SubStrType[indxTmp[!is.na(indxTmp)]], ' (subStructure.prob)')
    Comments(object) <- metIDcomments
    return(object)
}) # end function
