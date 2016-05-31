#' Identifies probable substructure type 
#' 
#' @description Identifies probable substructure type based on the summed
#' relative intensites and therefore proportion of total composite spectrum 
#' intensity explained.
#' @param a compMS2 object 
#' 
#' @return a data.frame of probable substructure annotations for each composite
#' spectrum, ranked by the sum of the relative intensities for that substructure
#' type 
#' @export
setGeneric("subStructure.prob", function(object, ...) standardGeneric("subStructure.prob"))

setMethod("subStructure.prob", signature = "CompMS2", function(object){
  
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if ("Frag.ID.type" %in% colnames(compSpectra(object)[[1]])){
    
    message("Identifying likely substructure type...")
    flush.console()
    
    names_CompSpec <- names(compSpectra(object))
    substrBestIds <- lapply(c(1:length(compSpectra(object))), function(x){
      spec.df <- compSpectra(object)[[x]]
      arrIndxAnno <- which(spec.df[, c("Frag.ID.type", "interfrag.loss.type", 
                                        "Neutral.loss.type")] != '', arr.ind = T)
      
      if(nrow(arrIndxAnno) > 0){
      subStrTypes <- spec.df[, c("Frag.ID.type", "interfrag.loss.type", 
                                 "Neutral.loss.type")]
      subStrTypes <- subStrTypes[arrIndxAnno]
      names(subStrTypes) <- spec.df$Rel_Intensity[arrIndxAnno[, 1]]
      relIntOrder <- order(spec.df$Rel_Intensity[arrIndxAnno[, 1]], decreasing = T)
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
      subStrs.tmp <- data.frame(table(subStrTypes), stringsAsFactors = F)
      subStrTypes <- tapply(as.numeric(names(subStrTypes)), subStrTypes, sum)
           # subStrTypes <- subStrTypes / sum(spec.df$intensity) 
      # summary results
      subStrTypes <- data.frame(compSpecName = rep(names_CompSpec[x], 
                                                   length(subStrTypes)),
                                SubStrType = names(subStrTypes), 
                                SumRelInt = subStrTypes,
                                stringsAsFactors = F)
      subStrTypes <- merge(subStrTypes, subStrs.tmp, by.x = "SubStrType", 
                           by.y = "subStrTypes")
      subStrTypes <- subStrTypes[order(subStrTypes$SumRelInt,
                                       decreasing = T), ]
      subStrTypes$nPeaks <- length(spec.df$intensity)
      } else {
      subStrTypes <- data.frame(compSpecName = names_CompSpec[x],
                                SubStrType = "no substructure detected", 
                                SumRelInt = "no substructure detected", 
                                Freq = "no substructure detected",
                                nPeaks = "no substructure detected",
                                stringsAsFactors = F)
      }
      return(subStrTypes)
    })
  } else {
    stop("subStructure.Annotate function has not yet been applied")
  }
     
    subStrAnno(object) <- do.call("rbind", substrBestIds)
    return(object)
    
  })
