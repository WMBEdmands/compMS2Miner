#'MS1featureMatch
#'
#' match ms1 features defined by mz and Rt to vector of ms2 features
#' 
#' @param mz mass to charge ratio (from the MS1 peak table).
#' @param RT retention time in seconds (from the MS1 peak table).
#' @param mz.v numeric vector of mass to charge ratio (from MS2 precursors)
#' @param RT.v  numeric vector of retention time in seconds (from MS2 precursors)
#' @param precursorPpm parts per million mass accuracy for match (default
#' @param ret retention time window (+/- seconds) for match
#' @param isoWid isolation width of ions.
#' @return returns indices of matches as a list object
#' 
MS1MatchSpectra <- function(metaData=NULL, MS2file=NULL, mz=NULL, RT=NULL, 
                            EIC=NULL, adduct=NULL, precursorPpm = 10, ret = 10, 
                            adducts=FALSE, isoWid=4){
  #error handling
  if(is.null(metaData)){
    stop("metaData value is missing with no default")
  } else if(is.null(MS2file)){
    stop("MS2file is missing with no default")
  } else if(is.null(mz)){
    stop("mz value is missing with no default")
  } else if(is.null(RT)){
    stop("RT value is missing with no default")
  } else if(is.null(EIC)){
    stop("EIC value is missing with no default")
  } else {
    # match MS1 mass by ppm tolerance to precursor m/z
    # match MS1 RT in seconds to precursor RT
    MS1.match <- which({metaData$precursorMz < 
                          mz+(mz/1E06)*precursorPpm &
                          metaData$precursorMz > 
                          mz-(mz/1E06)*precursorPpm &
                          metaData$retentionTime < RT+ret &
                          metaData$retentionTime > RT-ret &
                          metaData$MS2TICfilt.indx == 1} == TRUE)
    
    if(length(MS1.match) == 0){
      MS1.match <- 0
      return(MS1.match)
      
    } else  if(length(MS1.match) == 1){
      spectra <- data.frame(mzR::peaks(MS2file, MS1.match))
      colnames(spectra) <- c("mass", "intensity")
      #       spectra$scanSeqNum <- mzR::header(MS2file, MS1.match)$seqNum
      metaData.df <- metaData[MS1.match, , drop = FALSE]
      # ms1 scans
      precursorScans <- unique(metaData.df$precursorScanNum)
      # voltage ramps
      if(all(precursorScans == 0)){
        precursorScans <- metaData$acquisitionNum[MS1.match - 1]
        metaData.df$precursorScanNum <- precursorScans
      }
      ms1ScanIdx <- match(metaData$acquisitionNum, precursorScans)
      # subset if precursor scan missing
      metaData.df <- metaData.df[ms1ScanIdx[!is.na(ms1ScanIdx)], , drop=FALSE]
      precursorScans <- precursorScans[ms1ScanIdx[!is.na(ms1ScanIdx)]]
      ms1ScanIdx <- which(!is.na(ms1ScanIdx))
      ms1Prec <- unique(metaData.df$precursorMz)[1]
      isoChim <- data.frame(matrix(0, ncol=4, nrow=length(ms1ScanIdx)), stringsAsFactors = FALSE)
      colnames(isoChim) <-  c('isoMass', 'isoInt', 'possChim', 'maxInterIons')
      for(sc in 1:length(ms1ScanIdx)){
        x <- mzR::peaks(MS2file, ms1ScanIdx[sc])
        idxTmp <-  x[, 1] < {ms1Prec + 2.5} & x[, 1] > {ms1Prec - 2.5}
        if(sum(idxTmp) > 2){
        ms1SpecTmp <- x[idxTmp, , drop=FALSE]
        prIdx <- which.min(abs(ms1SpecTmp[, 1] - ms1Prec))
        iso1 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[prIdx, 1] + 1}))
        iso2 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[iso1, 1] + 1}))
        isoMass <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 1], 6), collapse = ';')
        isoInt <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 2], 2), collapse = ';')
        idxTmp <- x[, 1] < {ms1Prec + {isoWid/2}} & x[, 1] > {ms1Prec - {isoWid/2}} 
        possChim <- x[idxTmp, , drop=FALSE]
        prIdx <- which.min(abs(possChim[, 1] - ms1Prec))
        iso1 <- which.min(abs(possChim[, 1] - {possChim[prIdx, 1] + 1}))
        iso2 <- which.min(abs(possChim[, 1] - {possChim[iso1, 1] + 1}))
        possChim[, 2] <- {possChim[, 2]/possChim[prIdx, 2]} * 100
        # remove prec and isotopes within 0.5 Da
        idxClose <- possChim[, 1] < {possChim[prIdx, 1] + 0.5} & possChim[, 1] > {possChim[prIdx, 1] - 0.5} | possChim[, 1] < {possChim[iso1, 1] + 0.5} & possChim[, 1] > {possChim[iso1, 1] - 0.5} | possChim[, 1] < {possChim[iso2, 1] + 0.5} & possChim[, 1] > {possChim[iso2, 1] - 0.5}
        possChim <- possChim[idxClose == FALSE, , drop=FALSE]
        if(nrow(possChim) > 0){
          maxInterIons <- max(possChim[, 2])
          possChim <- any(possChim[, 2] >= 50)
        } else {
          maxInterIons <- 0
          possChim <- FALSE 
        }
        isoChim[sc, ] <- c(isoMass, isoInt, possChim, maxInterIons)
        }
      }
      row.names(isoChim) <- precursorScans
      metaData.df <- cbind(metaData.df, 
                           isoChim[match(metaData.df$precursorScanNum,
                           row.names(isoChim)), ])
      metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
      metaData.df$rtDiff <- metaData.df$retentionTime - RT
      metaData.df$MS1_EICno <- EIC
      metaData.df$MS1_mz <- mz
      metaData.df$MS1_RT <- RT
      
      if(adducts == TRUE){
      metaData.df$MS1_adduct <- adduct
      } else {
      metaData.df$MS1_adduct <- ''  
      }
      MS1.match <- list(spectra = spectra, metaData = metaData.df)    
    } else if(length(MS1.match) > 1){
      metaData.df <- metaData[MS1.match, , drop = FALSE]
      # ms1 scans
      precursorScans <- unique(metaData.df$precursorScanNum)
      if(all(precursorScans == 0)){
      precursorScans <- metaData$acquisitionNum[MS1.match - 1]
      metaData.df$precursorScanNum <- precursorScans
      } 
      ms1ScanIdx <- match(metaData$acquisitionNum, precursorScans)
      if(any(!is.na(ms1ScanIdx))){
      # subset if precursor scan missing
      metaData.df <- metaData.df[ms1ScanIdx[!is.na(ms1ScanIdx)], , drop=FALSE]
      precursorScans <- precursorScans[ms1ScanIdx[!is.na(ms1ScanIdx)]]
      MS1.match <- MS1.match[ms1ScanIdx[!is.na(ms1ScanIdx)]]
      if(length(MS1.match) == 1){
        spectra <- data.frame(mzR::peaks(MS2file, MS1.match))
      } else {
      spectra <- data.frame(do.call(rbind, mzR::peaks(MS2file, MS1.match)), stringsAsFactors = FALSE)
      }
      colnames(spectra) <- c("mass", "intensity")
      ms1ScanIdx <- which(!is.na(ms1ScanIdx))
      
      ms1Prec <- unique(metaData.df$precursorMz)[1]
      isoChim <- data.frame(matrix(NA, ncol=4, nrow=length(ms1ScanIdx)), stringsAsFactors = FALSE)
      colnames(isoChim) <-  c('isoMass', 'isoInt', 'possChim', 'maxInterIons')
      
      for(sc in 1:length(ms1ScanIdx)){
      x <- mzR::peaks(MS2file, ms1ScanIdx[sc])
      idxTmp <-  x[, 1] < {ms1Prec + 2.5} & x[, 1] > {ms1Prec - 2.5}
      if(sum(idxTmp) > 2){
      ms1SpecTmp <- x[idxTmp, , drop=FALSE]
      prIdx <- which.min(abs(ms1SpecTmp[, 1] - ms1Prec))
      iso1 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[prIdx, 1] + 1}))
      iso2 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[iso1, 1] + 1}))
      isoMass <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 1], 6), collapse = ';')
      isoInt <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 2], 2), collapse = ';')
      idxTmp <- x[, 1] < {ms1Prec + {isoWid/2}} & x[, 1] > {ms1Prec - {isoWid/2}} 
      possChim <- x[idxTmp, , drop=FALSE]
      prIdx <- which.min(abs(possChim[, 1] - ms1Prec))
      iso1 <- which.min(abs(possChim[, 1] - {possChim[prIdx, 1] + 1}))
      iso2 <- which.min(abs(possChim[, 1] - {possChim[iso1, 1] + 1}))
      possChim[, 2] <- {possChim[, 2]/possChim[prIdx, 2]} * 100
      # remove prec and isotopes
      idxClose <- possChim[, 1] < {possChim[prIdx, 1] + 0.5} & possChim[, 1] > {possChim[prIdx, 1] - 0.5} | possChim[, 1] < {possChim[iso1, 1] + 0.5} & possChim[, 1] > {possChim[iso1, 1] - 0.5} | possChim[, 1] < {possChim[iso2, 1] + 0.5} & possChim[, 1] > {possChim[iso2, 1] - 0.5}
      possChim <- possChim[idxClose == FALSE, , drop=FALSE]
      if(nrow(possChim) > 0){
      maxInterIons <- max(possChim[, 2])
      possChim <- any(possChim[, 2] >= 50)
      } else {
      maxInterIons <- 0
      possChim <- FALSE 
      }
      isoChim[sc, ] <- c(isoMass, isoInt, possChim, maxInterIons)
       }
      }
      row.names(isoChim) <- precursorScans
      metaData.df <- cbind(metaData.df, 
                           isoChim[match(metaData.df$precursorScanNum,
                                         row.names(isoChim)), ])     
      metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
      metaData.df$rtDiff <- metaData.df$retentionTime - RT
      metaData.df$MS1_EICno <- EIC
      metaData.df$MS1_mz <- mz
      metaData.df$MS1_RT <- RT
      
      if(adducts == TRUE){
        metaData.df$MS1_adduct <- adduct
      } else {
        metaData.df$MS1_adduct <- ''  
      }
      MS1.match <- list(spectra = spectra, metaData = metaData.df)    
      } else {
        MS1.match <- 0
      }
    }
    return(MS1.match)
  }
}
