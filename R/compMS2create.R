#' create a compMS2 object 
#'
#'@param MSfiles character vector of mzXML file locations
#'
#'
compMS2create <- function(MS2file = NULL, MS1features = NULL, 
                          TICfilter = 10000, precursorPpm = 10, ret = 10){
  
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
  # extract info from MS/MS file
  file.info.df <- MS2fileInfo(MS2file = MS2file, TICfilter = TICfilter)
  
  message("...DONE")
  flush.console()
  
  # cond if no MS2 level scans detected
  if(all(file.info.df$MS.scanType == 1)){
    # no MS2 level scans detected cond
    warning(paste0("No MS2 levels scans within ", MS2fileName, ",  check that the 
                   file has been converted to the mzXML format correctly."), 
            immediate.=T)
    flush.console()
    message("...moving to next MS2 file")
    flush.console()
  } else {
    # remove all MS/MS scans where the TIC is less than the minimum TIC threshold 
    # set by the user
    message(paste0("Of a total of ", length(which(file.info.df$MS.scanType == 2)),
                   " MS2 spectra..."))
    flush.console()
    # index ms2 scan and above TIC filter
    file.info.df$MS2TICfilt.indx <- (file.info.df$MS.scanType == 2 & 
                                     file.info.df$TICaboveFilter == 1) * 1
    nAboveTIC <- length(which(file.info.df$MS2TICfilt.indx == 1))
    message(paste0(nAboveTIC, " MS2 spectra were above the TIC filter of ", 
                   TICfilter))
    flush.console()
    # cond if no scan above the TIC filter
    if(length(nAboveTIC) == 0){ 
      warning(paste0("No MS2 levels scans above TIC filter of ", TICfilter, " in ", 
                     MS2fileName, ",  reduce the TIC filter parameter or check that 
                     the file has been converted to the mzXML format correctly."), 
              immediate.=T)
      flush.console()
      message("...moving to next MS2 file")
      flush.console()
    } else { 
      
      message("matching MS1 peak table to precursors of MS2 spectra...")
      flush.console()
      # mapply MS1 feature match
      MS1MS2match <- mapply(MS1MatchSpectra, EIC=MS1features[, 1], 
                            mz=MS1features[, 2], RT=MS1features[, 3], 
                            precursorPpm=precursorPpm, ret=ret, 
                            MoreArgs=list(file.info.df=file.info.df, 
                                            MS2file=MS2file))
      message("...done")
      flush.console()
      
      match.indx <- which(sapply(MS1MS2match, length) == 2)
      # calculate composite spectra
      message(paste0(length(match.indx), " MS1 features were matched to MS2 precursors"))
      flush.console()
      names(MS1MS2match) <- paste0(MS2fileName, "_", MS1features[, 1])
      MS1MS2match <- MS1MS2match[match.indx]
      return(MS1MS2match)
      #Results[names(MS1MS2match)] <- MS1MS2match
      close(MS2file)
      #time[i] <- (proc.time() - pmt)[["elapsed"]] 
    } # cond if no scan above the TIC filter
  } # cond if no MS2 level scans detected
} # end func