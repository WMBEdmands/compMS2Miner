#' combine ions across spectra matched to the same MS1 feature
#' @param mzError interpeak absolute m/z error for composite spectra signal grouping (default = 0.01)
#' @param minPeaks Minimum number of peaks per composite spectrum (default = 3)
#' @param specSimFilter umeric minimum spectral similarity score (dot product score) between spectra matched to the same MS1 feature (values between 0 to 1) if argument is supplied spectra are only combined if they have a minimum spectral similarity score (default = NULL). If all spectra matched to an MS1 feature are dissimilar from one another the MS2 spectrum with the highest precursor intensity will be returned.
#' @param binSizeMS2 numeric MS2 bin size for spectral similarity matching (default = 0.1) 
#' @export
setGeneric("combineMS2.Spectra", function(object, ...) standardGeneric("combineMS2.Spectra"))

setMethod("combineMS2.Spectra", signature = "CompMS2", function(object, 
                                                                mzError=0.01, 
                                                                minPeaks=3,
                                                                specSimFilter=NULL,
                                                                binSizeMS2=0.1){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  }
  # unique MS1 features vector
  MSfeatureNos <- gsub(".+_", "", names(object@compSpectra))
  
  if(!is.null(specSimFilter)){
  allSpecTmp <- do.call(rbind, lapply(object@compSpectra, function(x){
    massIntTmp <- x[, c('mass', 'intensity')]
   }))
  
 
  
  specNamesVecTmp <- unlist(mapply(rep, MSfeatureNos, each=sapply(object@compSpectra, nrow)))
  
  specNamesVecTmp <- paste0(specNamesVecTmp, '_', gsub('\\..+', '', row.names(allSpecTmp)))
  
  dotProdTmp <- dotProdMatrix(as.matrix(allSpecTmp), specNamesVecTmp, binSizeMS2=binSizeMS2)
  # calculate mean precursor intensity
  meanPrecInt <- sapply(object@metaData, function(x) mean(x$precursorIntensity))
  
 
  
  # spectra to combine
  specCombIndx <- unlist(lapply(unique(MSfeatureNos), function(x){
    grFeatTmp <- which(MSfeatureNos %in% x)
    dotPrTmp <- dotProdTmp[grFeatTmp, grFeatTmp, drop=F]
    indxTmp <- grFeatTmp[which(dotPrTmp[, 1] >= specSimFilter)]
    if(all(indxTmp == F)){
    indxTmp <- grFeatTmp[which.max(meanPrecInt[grFeatTmp])]  
    }
    # indxTmp[1] <- x
    return(indxTmp)
  }))

  message(nrow(dotProdTmp) - length(specCombIndx), ' spectra lower than the minimum spectral similarity score (', specSimFilter, ') were removed.\n')
  flush.console()
  
  specCombIndx <- sort(specCombIndx)
  object@compSpectra <- object@compSpectra[specCombIndx]
  object@metaData <- object@metaData[specCombIndx]
  MSfeatureNos <- MSfeatureNos[specCombIndx]
  } # end specSimFiltration
  
  ##############################################################################
    message(paste0("Combining ", length(object@compSpectra),
                   " spectra by MS1 feature number..."))
    flush.console()
    specGroups <-  split(object@compSpectra, as.factor(MSfeatureNos))
    # if specSimFilter supplied then calc spec Sim
    specGroups <- lapply(specGroups, function(x) do.call(rbind, x))
    metaDataGroups <-  tapply(object@metaData, as.factor(MSfeatureNos), function(x) rbind(x))
    metaDataGroups <- lapply(metaDataGroups, function(x){ 
      tmp <- do.call(c, x)
      tmp.names <-  paste(rep(colnames(x), 
                              each = length(tmp)/ ncol(x)),
                          names(tmp), sep = "_")
      names(tmp) <- tmp.names
      return(tmp)
    })
    
    if(object@Parameters$nCores > 0){
      
      # create a cluster using the doSNOW package
      message(paste0("Starting SNOW cluster with ", object@Parameters$nCores,
                     " local sockets..."))
      flush.console()
      
      cl <- parallel::makeCluster(object@Parameters$nCores) 
      doSNOW::registerDoSNOW(cl)
      
      # foreach and dopar from foreach package
      sign.group <- foreach(j = 1:length(specGroups),
                            .packages = c('stats')) %dopar% {
                              signalGrouping(spectrum.df = specGroups[[j]], 
                                             mzError=mzError, 
                                             minPeaks = minPeaks)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
    } else {
      # create list to store results
      sign.group <- vector("list", length(specGroups))
      # create progress bar
      pb <- txtProgressBar(min=0, max=length(sign.group), style=3)
      
      for(j in 1:length(sign.group)){
        
        #progress bar
        Sys.sleep(0.01)
        setTxtProgressBar(pb, j)
        flush.console()
        
        sign.group.tmp <- signalGrouping(spectrum.df = specGroups[[j]], 
                                         mzError=mzError, 
                                         minPeaks = minPeaks)
        
        sign.group[[j]] <- sign.group.tmp
      }
    }
    
    # logical if no peaks returned
    groupIndx <- sapply(sign.group, function(x) !is.character(x))
    
    message("...done")
    flush.console()
    # number of comp spectra returned
    message(sum(groupIndx), " composite spectra contained more than or equal to ",
            minPeaks," peaks following ion grouping")
    flush.console()
    # return grouped
    # calculate number of interfragment difference lower than 0.1 m/z and inform
    # user
    nInterFragless <- sapply(sign.group[groupIndx], function(x){
      intfrag.diff <- as.numeric(c(diff(x[, "mass"]), 0))
      length(which(intfrag.diff < 0.1))})
    message("The range of interfragment differences less than 0.1 m/z in the composite spectra is ",
            paste(c("min :", " max :"), range(nInterFragless)), "\n")
    flush.console()
    message("The average number of interfragment differences less than 0.1 m/z in the composite spectra is ",
            round(mean(nInterFragless), digits=0), "\n")
    flush.console()
    
    if(round(mean(nInterFragless), digits=0) > 2){
      warning("The average number of interfragment differences less than 0.1 m/z is greater than 2: please consider increasing the mzError parameter above ", 
              mzError)
    }
    
    names.tmp <- paste0("CC_",names(specGroups))
    names(sign.group) <- names.tmp
    compSpectra(object) <- sign.group[groupIndx]
    names(metaDataGroups) <- names.tmp
    metaData(object) <- metaDataGroups[groupIndx]
    
    return(object)
  
}) # end function
