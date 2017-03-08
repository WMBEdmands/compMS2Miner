#' create a compMS2 object 
#'
#'@param MSfiles character vector of mzXML file locations
compMS2Create <- function(MS2file = NULL, MS1features = NULL, 
                          TICfilter = 10000, precursorPpm = 10, ret = 10, 
                          adducts=FALSE, isoWid=4){
  
  # MS2 file name
  MS2fileName <- basename(MS2file)
  
  message(paste0("Reading ", MS2fileName, "..."))
  flush.console()
  # read MS2 file
  MS2file <- mzR::openMSfile(MS2file)
  
  message("...DONE")
  flush.console()
  
  message("extracting metaData from MS2 file")
  flush.console()
  
  metaData <- mzR::header(MS2file)
  metaData <- metaData[, c('msLevel', 'precursorMZ', 'retentionTime',
                           'totIonCurrent', 'precursorIntensity',
                           'collisionEnergy', 'basePeakMZ', 'basePeakIntensity',
                           'acquisitionNum', 'precursorScanNum')]
  metaData$TICaboveFilter <- {metaData$totIonCurrent >= TICfilter} * 1
  colnames(metaData) <- c("MS.scanType", "precursorMz", "retentionTime", "TIC", 
                          "precursorIntensity", 
                          "collisionEnergy", "basePeakMz", "basePeakIntensity",
                          "acquisitionNum","precursorScanNum", "TICaboveFilter")
  
  metaData <- metaData[, c("MS.scanType", "precursorMz", "retentionTime", "TIC", 
                           "TICaboveFilter", "precursorIntensity", 
                           "collisionEnergy", "basePeakMz", "basePeakIntensity",
                           "acquisitionNum","precursorScanNum")]
  message("...DONE")
  flush.console()
  
  # cond if no MS2 level scans detected
  if(all(metaData$MS.scanType == 1)){
    # no MS2 level scans detected cond
    warning(paste0("No MS2 levels scans within ", MS2fileName, ",  check that the 
                   file has been converted to the mzXML format correctly."), 
            immediate.=TRUE)
    flush.console()
    message("...moving to next MS2 file")
    flush.console()
  } else {
    # remove all MS/MS scans where the TIC is less than the minimum TIC threshold 
    # set by the user
    message(paste0("Of a total of ", length(which(metaData$MS.scanType == 2)),
                   " MS2 spectra..."))
    flush.console()
    # index ms2 scan and above TIC filter
    metaData$MS2TICfilt.indx <- (metaData$MS.scanType == 2 & 
                                     metaData$TICaboveFilter == 1) * 1
    nAboveTIC <- length(which(metaData$MS2TICfilt.indx == 1))
    message(paste0(nAboveTIC, " MS2 spectra were above the TIC filter of ", 
                   TICfilter))
    flush.console()
    # cond if no scan above the TIC filter
    if(length(nAboveTIC) == 0){ 
      warning(paste0("No MS2 levels scans above TIC filter of ", TICfilter, " in ", 
                     MS2fileName, ",  reduce the TIC filter parameter or check that 
                     the file has been converted to the mzXML format correctly."), 
              immediate.=TRUE)
      flush.console()
      message("...moving to next MS2 file")
      flush.console()
    } else { 
      
      message("matching MS1 peak table to precursors of MS2 spectra...")
      flush.console()
      # mapply MS1 feature match
      MS1MS2match <- mapply(MS1MatchSpectra, EIC=MS1features[, 1], 
                            mz=MS1features[, 2], RT=MS1features[, 3], 
                            adduct=MS1features[, 4],
                            precursorPpm=precursorPpm, ret=ret, 
                            MoreArgs=list(metaData=metaData, 
                                          MS2file=MS2file, adducts=adducts, 
                                          isoWid=isoWid))
      
      # for(i in 1:nrow(MS1features)){
      #   tmp <- MS1MatchSpectra(EIC=MS1features[i, 1], 
      #                          mz=MS1features[i, 2], RT=MS1features[i, 3], 
      #                          adduct=MS1features[i, 4],
      #                          precursorPpm=precursorPpm, ret=ret,
      #                          metaData=metaData, 
      #                          MS2file=MS2file, adducts=adducts, 
      #                          isoWid=isoWid)
      #   # EIC=MS1features[i, 1]; 
      #   # mz=MS1features[i, 2]; RT=MS1features[i, 3]; 
      #   # adduct=MS1features[i, 4];
      # }
      message("...done")
      flush.console()
      
      match.indx <- which(sapply(MS1MS2match, length) == 2)
      # calculate composite spectra
      message(paste0(length(match.indx), " MS1 features were matched to MS2 precursors"))
      flush.console()
      names(MS1MS2match) <- paste0(MS2fileName, "_", MS1features[, 1])
      
      # check for chimeric spectra and isotopes
      
      
      MS1MS2match <- MS1MS2match[match.indx]
      return(MS1MS2match)
      #Results[names(MS1MS2match)] <- MS1MS2match
      # close(MS2file)
      #time[i] <- (proc.time() - pmt)[["elapsed"]] 
    } # cond if no scan above the TIC filter
  } # cond if no MS2 level scans detected
} # end func
