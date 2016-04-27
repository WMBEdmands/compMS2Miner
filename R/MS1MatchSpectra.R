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
#' 
#' @return returns indices of matches as a list object
#' 
MS1MatchSpectra <- function(file.info.df=NULL, MS2file=NULL, mz=NULL, RT=NULL, 
                            EIC=NULL, precursorPpm = 10, ret = 10){
  #error handling
  if(is.null(file.info.df)){
    stop("file.info.df value is missing with no default")
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
    MS1.match <- which((file.info.df$precursorMz < 
                          mz+(mz/1E06)*precursorPpm &
                          file.info.df$precursorMz > 
                          mz-(mz/1E06)*precursorPpm &
                          file.info.df$retentionTime < RT+ret &
                          file.info.df$retentionTime > RT-ret &
                          file.info.df$MS2TICfilt.indx == 1) == T)
    
    if(length(MS1.match) == 0){
      MS1.match <- 0
      return(MS1.match)
      
    } else  if(length(MS1.match) == 1){
      spectra <- data.frame(mzR::peaks(MS2file, MS1.match))
      colnames(spectra) <- c("mass", "intensity")
      #       spectra$scanSeqNum <- mzR::header(MS2file, MS1.match)$seqNum
      metaData.df <- file.info.df[MS1.match, , drop = F]
      #       spectra$precursorScanNum <- metaData.df$precursorScanNum
      metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
      metaData.df$rtDiff <- metaData.df$retentionTime - RT
      metaData.df$MS1_EICno <- EIC
      metaData.df$MS1_mz <- mz
      metaData.df$MS1_RT <- RT
      MS1.match <- list(spectra = spectra, metaData = metaData.df)    
    } else if(length(MS1.match) > 1){
      spectra <- data.frame(do.call(rbind, mzR::peaks(MS2file, MS1.match)))
      colnames(spectra) <- c("mass", "intensity")
      metaData.df <- file.info.df[MS1.match, , drop = F]
      metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
      metaData.df$rtDiff <- metaData.df$retentionTime - RT
      metaData.df$MS1_EICno <- EIC
      metaData.df$MS1_mz <- mz
      metaData.df$MS1_RT <- RT
      MS1.match <- list(spectra = spectra, metaData = metaData.df)    
    }
    return(MS1.match)
  }
}