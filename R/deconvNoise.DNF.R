#' Spectral noise filtration using dynamic noise filter
#' 
#' @description uses the dynamic noise filtration algorithm adapted from the method described 
#' in Xu H. and Frietas M. "A Dynamic Noise Level Algorithm for Spectral 
#' Screening of Peptide MS/MS Spectra" 2010 BMC Bioinformatics. 

#' @param DNF. numeric dynamic noise filter minimum signal to noise threshold 
#' (default = 2), calculated as the ratio between the linear model predicted 
#' intensity value and the actual intensity.
#' @param minPeaks. integer minimum number of signal peaks following dynamic 
#' noise filtration (default = 1).
#' @param maxPeaks. integer maximum number of signal peaks the function will continue
#' until both the minimum DNF signal to noise ratio is exceeding and the number
#' of peaks is lower than the maximum (default = 60).
#' @param minInt. numeric minimum intensity to commence the dynamic noise filter
#' algorithm. Low values will increase computation time and increase the chance
#' that the DNF algorithm will terminate prematurely (default = 250).
#' @return noise filtered MS2 spectra.
#' @param verbose logical if TRUE display progress bars.
#' @source Xu H. and Frietas M. "A Dynamic Noise Level Algorithm for Spectral 
#' Screening of Peptide MS/MS Spectra" 2010 BMC Bioinformatics.
#' @export 
setGeneric("deconvNoise.DNF", function(object, ...) standardGeneric("deconvNoise.DNF"))

setMethod("deconvNoise.DNF", signature = "compMS2", function(object, DNF=2, 
                                                             minPeaks=1,
                                                             maxPeaks=20,
                                                             minInt=250, 
                                                             verbose=TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    message(paste0("Applying dynamic noise filter to ", length(compSpectra(object)),
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
      noise.filt <- foreach(j = 1:length(compSpectra(object)),
                            .packages = c('RcppEigen', 'compMS2Miner'), .options.snow=opts) %dopar% {
#                               for(j in 1:length(compSpectra(object))){
#                               compMS2Miner::dynamicNoiseFilter(spectrum.df=compSpectra(object)[[j]], 
                              dynamicNoiseFilter(spectrum.df=compSpectra(object)[[j]], 
                                                 DNF=DNF, minPeaks=minPeaks, 
                                                 maxPeaks=maxPeaks, 
                                                 minInt=minInt)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
    } else {
      # create list to store results
      noise.filt <- vector("list", length(compSpectra(object)))
      # create progress bar
      if(verbose == TRUE){ pb <- txtProgressBar(min=0, max=length(noise.filt), style=3)}
      
      for(j in 1:length(noise.filt)){
        
        #progress bar
        if(verbose==TRUE){setTxtProgressBar(pb, j)}
        flush.console()
        
        noise.filt.tmp <- dynamicNoiseFilter(spectrum.df=compSpectra(object)[[j]], 
                                             DNF=DNF, minPeaks=minPeaks, 
                                             maxPeaks=maxPeaks, minInt=minInt)
        
        noise.filt[[j]] <- noise.filt.tmp
      }
    }
    # logical if no peaks returned
    noiseIndx <- sapply(noise.filt, function(x) x$aboveMinPeak == TRUE)
    
    message("...done")
    flush.console()
    # number of comp spectra returned
    message(sum(noiseIndx), " spectra contained more than or equal to ",
            minPeaks," peaks following dynamic noise filtration")
    flush.console()
    # return noise filtered
    compSpec.tmp <- lapply(noise.filt, function(x) x$Above.noise)
    prevMetaData.tmp <- metaData(object)
    metaData.tmp <- lapply(c(1:length(noise.filt)), function(x) 
      c(prevMetaData.tmp[[x]],
        data.frame(noise.filt[[x]]$metaData, stringsAsFactors = FALSE)))
    
    names.tmp <- names(compSpectra(object))
    names(compSpec.tmp) <- names.tmp
    names(metaData.tmp) <- names.tmp
    compSpectra(object) <- compSpec.tmp[noiseIndx]
    metaData(object) <- metaData.tmp[noiseIndx]
    
    return(object)
    
  }
}) # end function

