#' Spectral noise filtration using dynamic noise filter
#' 
#' @description uses the dynamic noise filtration algorithm adapted from the method described 
#' in Xu H. and Frietas M. "A Dynamic Noise Level Algorithm for Spectral 
#' Screening of Peptide MS/MS Spectra" 2010 BMC Bioinformatics. 

#' @param DNF. dynamic noise filter minimum signal to noise threshold 
#' (default = 2), calculated as the ratio between the linear model predicted 
#' intensity value and the actual intensity.
#' @param minPeaks. minimum number of signal peaks following dynamic 
#' noise filtration (default = 5).
#' @param maxPeaks. maximum number of signal peaks the function will continue
#' until both the minimum DNF signal to noise ratio is exceeding and the number
#' of peaks is lower than the maximum (default = 5).
#' 
#' @return noise filtered composite MS2 spectra.
#' @source Xu H. and Frietas M. "A Dynamic Noise Level Algorithm for Spectral 
#' Screening of Peptide MS/MS Spectra" 2010 BMC Bioinformatics.
#' @export 
setGeneric("deconvNoise.DNF", function(object, ...) standardGeneric("deconvNoise.DNF"))

setMethod("deconvNoise.DNF", signature = "CompMS2", function(object, DNF=2, 
                                                             minPeaks=5,
                                                             maxPeaks=20,
                                                             minInt=250){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    message(paste0("Applying dynamic noise filter to ", length(object@compSpectra),
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
      noise.filt <- foreach(j = 1:length(object@compSpectra),
                            .packages = c('RcppEigen')) %dopar% {
#                               for(j in 1:length(object@compSpectra)){
#                               CompMS2miner::dynamicNoiseFilter(spectrum.df=object@compSpectra[[j]], 
                              dynamicNoiseFilter(spectrum.df=object@compSpectra[[j]], 
                                                 DNF=DNF, minPeaks=minPeaks, 
                                                 maxPeaks=maxPeaks, 
                                                 minInt=minInt)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
    } else {
      # create list to store results
      noise.filt <- vector("list", length(object@compSpectra))
      # create progress bar
      pb <- txtProgressBar(min=0, max=length(noise.filt), style=3)
      
      for(j in 1:length(noise.filt)){
        
        #progress bar
        Sys.sleep(0.01)
        setTxtProgressBar(pb, j)
        flush.console()
        
        noise.filt.tmp <- dynamicNoiseFilter(spectrum.df=object@compSpectra[[j]], 
                                             DNF=DNF, minPeaks=minPeaks, 
                                             maxPeaks=maxPeaks, minInt=minInt)
        
        noise.filt[[j]] <- noise.filt.tmp
      }
    }
    # logical if no peaks returned
    noiseIndx <- sapply(noise.filt, function(x) x$aboveMinPeak == T)
    
    message("...done")
    flush.console()
    # number of comp spectra returned
    message(sum(noiseIndx), " composite spectra contained more than or equal to ",
            minPeaks," peaks following dynamic noise filtration")
    flush.console()
    # return noise filtered
    compSpec.tmp <- lapply(noise.filt, function(x) x$Above.noise)
    prevMetaData.tmp <- metaData(object)
    metaData.tmp <- lapply(c(1:length(noise.filt)), function(x) 
      c(prevMetaData.tmp[[x]],
        data.frame(noise.filt[[x]]$metaData)))
    
    names.tmp <- names(compSpectra(object))
    names(compSpec.tmp) <- names.tmp
    names(metaData.tmp) <- names.tmp
    compSpectra(object) <- compSpec.tmp[noiseIndx]
    metaData(object) <- metaData.tmp[noiseIndx]
    
    return(object)
    
  }
})