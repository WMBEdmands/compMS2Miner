#' combine ions across spectra matched to the same MS1 feature
#' @param mzError interpeak absolute m/z error for composite spectra signal grouping (default = 0.01)
#' @param minPeaks Minimum number of peaks per composite spectrum (default = 1)
#' @param specSimFilter numeric minimum spectral similarity score (dot product score) between spectra matched to the same MS1 feature (values between 0 to 1) if argument is supplied spectra are only combined if they have a minimum spectral similarity score (default = NULL). If all spectra matched to an MS1 feature are dissimilar from one another the MS2 spectrum with the highest precursor intensity will be returned.
#' @param binSizeMS2 numeric MS2 bin size for spectral similarity matching (default = 0.1) 
#' @param verbose logical if TRUE display progress bars.
#' @export
setGeneric("combineMS2.Spectra", function(object, ...) standardGeneric("combineMS2.Spectra"))

setMethod("combineMS2.Spectra", signature = "compMS2", function(object, 
                                                                mzError=0.01, 
                                                                minPeaks=1,
                                                                specSimFilter=NULL,
                                                                binSizeMS2=0.1,
                                                                verbose=TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  }
  # unique MS1 features vector
  MSfeatureNos <- gsub(".+_", "", names(compSpectra(object)))
  ##############################################################################
    message(paste0("Combining ", length(compSpectra(object)),
                   " spectra by MS1 feature number..."))
    flush.console()
    specGroups <-  split(compSpectra(object), as.factor(MSfeatureNos))
    # if specSimFilter supplied then calc spec Sim
    specGroups <- lapply(specGroups, function(x){
                         specTmp <- do.call(rbind, x)
                         specTmp <- cbind(specTmp, gsub('\\.[0-9]+$', '', row.names(specTmp)))})
   
    
    if(Parameters(object)$nCores > 0){
      
      # create a cluster using the doSNOW package
      message(paste0("Starting SNOW cluster with ", Parameters(object)$nCores,
                     " local sockets..."))
      flush.console()
      
      cl <- parallel::makeCluster(Parameters(object)$nCores, outfile='') 
      doSNOW::registerDoSNOW(cl)
      
      progSeq <- round({length(specGroups) * seq(0, 1, 0.05)}, 0)
      progSeq[1] <- 1
      cat(paste0('Progress (', length(specGroups), ' spectrum groups):\n'))
      progress <- function(n){if(n %in% progSeq){cat(paste0(round({n/length(specGroups)} * 100, 0), '%  '))}}
      if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
      
      # foreach and dopar from foreach package
      sign.group <- foreach(j = 1:length(specGroups),
                            .packages = c('stats'), .options.snow=opts) %dopar% {
                              spectrum.df <- specGroups[[j]]
                              if(!is.null(specSimFilter)){
                                dotProdTmp <- suppressMessages(dotProdMatrix(as.matrix(spectrum.df[, 1:2]), spectrum.df[, 3], binSizeMS2=binSizeMS2)) 
                                indxTmp <- which(dotProdTmp[, 1] >= specSimFilter)
                                if(length(indxTmp) == 0){
                                  # if neither similar take the spectrum with the highest TIC
                                  indxTmp <- which.max(tapply(spectrum.df[, 2], spectrum.df[, 3], sum))  
                                }
                                # subset spectrum
                                retSpecIndx <- unique(spectrum.df[, 3])[indxTmp]
                                retSpecIndx <- spectrum.df[, 3] %in% retSpecIndx
                                spectrum.df <- spectrum.df[retSpecIndx, , drop=FALSE]
                              }
                              remSpec <- unique(spectrum.df[, 3])
                              spectrum.df <- spectrum.df[, 1:2]
                              list(spec=signalGrouping(spectrum.df = spectrum.df, 
                                             mzError=mzError, 
                                             minPeaks = minPeaks), remSpec=remSpec)}
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
    } else {
      # create list to store results
      sign.group <- vector("list", length(specGroups))
      # create progress bar
      if(verbose == TRUE){ pb <- txtProgressBar(min=0, max=length(sign.group), style=3)}
      
      for(j in 1:length(sign.group)){
        
        #progress bar
        Sys.sleep(0.01)
        if(verbose==TRUE){setTxtProgressBar(pb, j)}
        flush.console()
        spectrum.df <- specGroups[[j]]
        if(!is.null(specSimFilter)){
          dotProdTmp <- suppressMessages(dotProdMatrix(as.matrix(spectrum.df[, 1:2, drop=FALSE]), spectrum.df[, 3], binSizeMS2=binSizeMS2)) 
          indxTmp <- which(dotProdTmp[, 1] >= specSimFilter)
          if(length(indxTmp) == 0){
            # if neither similar take the spectrum with the highest TIC
            indxTmp <- which.max(tapply(spectrum.df[, 2], spectrum.df[, 3], sum))  
          }
          # subset spectrum
          retSpecIndx <- unique(spectrum.df[, 3])[indxTmp]
          retSpecIndx <- spectrum.df[, 3] %in% retSpecIndx
          spectrum.df <- spectrum.df[retSpecIndx, , drop=FALSE]
        }
        remSpec <- unique(spectrum.df[, 3])
        spectrum.df <- spectrum.df[, 1:2]
        
        sign.group.tmp <- list(spec=signalGrouping(spectrum.df = spectrum.df, 
                                         mzError=mzError, 
                                         minPeaks = minPeaks), remSpec=remSpec)
        
        sign.group[[j]] <- sign.group.tmp
      }
    }
    # remaining spectra
    remainSpec <- lapply(sign.group, function(x){x$remSpec})
    remainSpec <- unlist(remainSpec)
    if(length(remainSpec) != length(compSpectra(object))){
      message(length(compSpectra(object)) - length(remainSpec), ' spectra lower than the minimum spectral similarity score (', specSimFilter, ') were removed.\n')
      flush.console()
      indxTmp <- names(metaData(object)) %in% remainSpec
      metaData(object) <- metaData(object)[indxTmp]
      compSpectra(object) <- compSpectra(object)[indxTmp]
      MSfeatureNos <- MSfeatureNos[indxTmp]
    }
    sign.group <- lapply(sign.group, function(x) x[[1]])
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
    
    metaDataGroups <- vector('list', length(metaData(object)))
    for(k in 1:length(metaData(object))){
      namesTmp <- paste0(names(metaData(object)[k]), '_', names(metaData(object)[[k]]))
      metaDataGroups[[k]] <- metaData(object)[[k]]
      names(metaDataGroups[[k]]) <- namesTmp
    }
    
    metaDataGroups <- lapply(split(metaDataGroups, as.factor(MSfeatureNos)), function(x) do.call(c, x))
    names.tmp <- paste0("CC_",names(specGroups))
    names(sign.group) <- names.tmp
    compSpectra(object) <- sign.group[groupIndx]
    names(metaDataGroups) <- names.tmp
    metaData(object) <- metaDataGroups[groupIndx]
    
    return(object)
  
}) # end function
