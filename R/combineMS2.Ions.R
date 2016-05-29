#' Combine spectra peaks within a composite spectra
#' 
#' @details group ions according to absolute m/z error. 
#' The default parameters are suitable for a high-resolution Q-ToF.
#' Following ion grouping, signal intensities are summed and an average m/z 
#' calculated for each ion group. This signal summing serves to increase the
#' overall intensity of true ion signal across multiple scans and reduce the 
#' contribution of noise within the spectrum. Calculation of the central tendency
#' of each ion group serve to homogenize the random error and improve the mass
#' accuracy of each composite spectrum peak.    
#' 
#' @param mzError interpeak absolute m/z error for composite spectra signal grouping.
#' @param minPeaks Minimum number of peaks per composite spectrum. 
#' @param ... option arguments to be passed along.
#' 
#' @return A compMS2 object with ion grouped composite spectra.
#' @export
setGeneric("combineMS2.Ions", function(object, ...) standardGeneric("combineMS2.Ions"))

setMethod("combineMS2.Ions", signature = "CompMS2", function(object, 
                                                             mzError=0.001, 
                                                             minPeaks=5){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    message(paste0("Grouping ions in ", length(object@compSpectra),
                   " composite spectra..."))
    flush.console()
    
    if(object@Parameters$nCores > 0){
      # create a cluster using the doSNOW package
      message(paste0("Starting SNOW cluster with ", object@Parameters$nCores,
                     " local sockets..."))
      flush.console()
      
      cl <- parallel::makeCluster(object@Parameters$nCores) 
      doSNOW::registerDoSNOW(cl)
      
      # foreach and dopar from foreach package
      sign.group <- foreach(j = 1:length(object@compSpectra),
                            .packages = c('stats')) %dopar% {
                              signalGrouping(spectrum.df = object@compSpectra[[j]], 
                                             mzError=mzError,
                                             minPeaks=minPeaks)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
      
    } else {
      # create list to store results
      sign.group <- vector("list", length(object@compSpectra))
      # create progress bar
      pb <- txtProgressBar(min=0, max=length(sign.group), style=3)
      
      for(j in 1:length(sign.group)){
        
        #progress bar
        Sys.sleep(0.01)
        setTxtProgressBar(pb, j)
        flush.console()
        
        sign.group.tmp <- signalGrouping(spectrum.df = object@compSpectra[[j]], 
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
    
    # return grouped
    
    names.tmp <- names(compSpectra(object))
    names(sign.group) <- names.tmp
    compSpectra(object) <- sign.group[groupIndx]
    metaData(object) <- metaData(object)[groupIndx]
    
    return(object)
  }
})