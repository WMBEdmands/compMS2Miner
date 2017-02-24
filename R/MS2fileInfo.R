#' MS2 file information
#' 
#' Extract precursor mass-to-charge ratio, retention time, scan type and total
#' ion current for each MS2 file scan.
#' 
#' @param MS2file MS2 file imported into R via the readMzXml package as a list 
#' @return data frame number of rows equal to the number of scans and 6 
#' observations:
#' 1. "MS.scanType" = MS scan type, 1 or 2.
#' 2. "precursorMz" = numeric mass-to-charge ratio for each MS2 scan
#' 3. "retentionTime" = numeric precusor retention time for each MS2 scan
#' 4. "TIC" = Total ion current for this scan.
#' 5. "TICaboveFilter" = boolean if TIC above minimum TIC filter equal 1, else 0
#' 6. "precursorIntensity" = MS2 precursor intensity, if MS1 scan returns
#' zero
#' 7. "collisionEnergy" = collision energy (eV)
#' 8. "basePeakMz" = mass-to-charge ratio of the base peak for the scan
#' 9. "basePeakIntensity" = intensity of the base peak for the scan 
#' 10. "precursorScanNum" = precursor scan number (MS2) or scan number (MS1) 
#' @export
MS2fileInfo <- function(MS2file=NULL, TICfilter=NULL){
  # error handling 
  if(is.null(TICfilter)){
    
    stop("A TICfilter value must be supplied (e.g. 10,000)")
    
  } else if(is.null(MS2file)){
    
    stop("an MS2file (.mzXML) object must be supplied")
    # error handling to if mzXML file properly formed
  } else {
    
    Info<-function(x){
      tmp.info <- c((mzR::header(MS2file, x)$msLevel == 2)+1,
                    round(mzR::header(MS2file, x)$precursorMZ, digits = 6),
                    round(mzR::header(MS2file, x)$retentionTime, digits = 4),  
                    mzR::header(MS2file, x)$totIonCurrent,
                    (mzR::header(MS2file, x)$totIonCurrent >= TICfilter)*1,
                    mzR::header(MS2file, x)$precursorIntensity,
                    round(mzR::header(MS2file, x)$collisionEnergy, digits = 1),
                    round(mzR::header(MS2file, x)$basePeakMZ, digits=4),
                    round(mzR::header(MS2file, x)$basePeakIntensity, digits=1),
                    mzR::header(MS2file, x)$acquisitionNum, 
                    mzR::header(MS2file, x)$precursorScanNum)}    
    
    # create data frame using info function     
    tmp.df<-as.data.frame(t(sapply(c(1:length(MS2file)),FUN=Info)),
                          stringsAsFactors = FALSE) 
    # add column names
    colnames(tmp.df) <- c("MS.scanType", "precursorMz", "retentionTime", "TIC", 
                          "TICaboveFilter", "precursorIntensity", 
                          "collisionEnergy", "basePeakMz", "basePeakIntensity",
                          "acquisitionNum","precursorScanNum")
    
    return(tmp.df)
  }
}
