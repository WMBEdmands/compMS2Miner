#' Combine spectra peaks within individual spectra
#' 
#' @details group ions according to absolute m/z error. 
#' The default parameters are suitable for a high-resolution Q-ToF.
#' Following ion grouping, signal intensities are summed and an average m/z 
#' calculated for each ion group. This signal summing serves to increase the
#' overall intensity of true ion signal across multiple scans and reduce the 
#' contribution of noise within the spectrum. Calculation of the central tendency
#' of each ion group serve to homogenize the random error and improve the mass
#' accuracy of each spectrum peak.    
#' 
#' @param mzError interpeak absolute m/z error for spectra signal grouping (default = 0.01).
#' @param minPeaks Minimum number of peaks per spectrum (default = 1). 
#' @param ... option arguments to be passed along.
#' @param verbose logical if TRUE display progress bars.
#' @return A compMS2 object with ion grouped composite spectra.
#' @export
setGeneric("combineMS2.Ions", function(object, ...) standardGeneric("combineMS2.Ions"))

setMethod("combineMS2.Ions", signature = "compMS2", function(object, 
                                                             mzError=0.01, 
                                                             minPeaks=1,
                                                             verbose=TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    message(paste0("Grouping ions in ", length(compSpectra(object)),
                   " spectra..."))
    flush.console()
    
    if(Parameters(object)$nCores > 0){
      if(!require(foreach)){
        stop('package foreach must be installed to use this function in parallel')
      }
      if(!require(doSNOW)){
        stop('package doSNOW must be installed to use this function in parallel')
      }
      # create a cluster using the doSNOW package
      message(paste0("Starting SNOW cluster with ", Parameters(object)$nCores,
                     " local sockets..."))
      flush.console()
      
      cl <- parallel::makeCluster(Parameters(object)$nCores, outfile='') 
      doSNOW::registerDoSNOW(cl)
      progSeq <- round({length(compSpectra(object)) * seq(0, 1, 0.05)}, 0)
      progSeq[1] <- 1
      cat(paste0('Progress (', length(compSpectra(object)), ' spectra):\n'))
      progress <- function(n){if(n %in% progSeq){cat(paste0(round({n/length(compSpectra(object))} * 100, 0), '%  '))}}
      if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
      # foreach and dopar from foreach package
      sign.group <- foreach(j = 1:length(compSpectra(object)),
                            .packages = c('stats'), .options.snow=opts) %dopar% {
                              signalGrouping(spectrum.df = compSpectra(object)[[j]], 
                                             mzError=mzError,
                                             minPeaks=minPeaks)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
      
    } else {
      # create list to store results
      sign.group <- vector("list", length(compSpectra(object)))
      # create progress bar
      if(verbose == TRUE){ pb <- txtProgressBar(min=0, max=length(sign.group), style=3)}
      
      for(j in 1:length(sign.group)){
        
        #progress bar
        if(verbose==TRUE){setTxtProgressBar(pb, j)}
        flush.console()
        
        sign.group.tmp <- signalGrouping(spectrum.df = compSpectra(object)[[j]], 
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
    message(sum(groupIndx), " spectra contained more than or equal to ",
            minPeaks," peaks following ion grouping")
    flush.console()
    
    # calculate number of interfragment difference lower than 0.1 m/z and inform
    # user
    nInterFragless <- sapply(sign.group[groupIndx], function(x){
      intfrag.diff <- as.numeric(c(diff(x[, "mass"]), 0))
      length(which(intfrag.diff < 0.1))})
      message("The range of interfragment differences less than 0.1 m/z in the spectra is ",
              paste(c("min :", " max :"), range(nInterFragless)), "\n")
      flush.console()
      message("The average number of interfragment differences less than 0.1 m/z in the spectra is ",
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
