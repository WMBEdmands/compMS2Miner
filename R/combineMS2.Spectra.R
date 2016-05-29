#' combine ions across spectra matched to the same MS1 feature
#' @param mzError interpeak absolute m/z error for composite spectra signal grouping
#' @param minPeaks Minimum number of peaks per composite spectrum 
#' @export
setGeneric("combineMS2.Spectra", function(object, ...) standardGeneric("combineMS2.Spectra"))

setMethod("combineMS2.Spectra", signature = "CompMS2", function(object, 
                                                                mzError=0.001, 
                                                                minPeaks = 5){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    
    message(paste0("Combining ", length(object@compSpectra),
                   " composite spectra by MS1 feature number..."))
    flush.console()
    MSfeature.nums <- gsub(".+_", "", names(object@compSpectra))
    specGroups <-  tapply(object@compSpectra, as.factor(MSfeature.nums), function(x) rbind(x))
    specGroups <- lapply(specGroups, function(x) do.call(rbind, x))
    metaDataGroups <-  tapply(object@metaData, as.factor(MSfeature.nums), function(x) rbind(x))
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
  }
})
