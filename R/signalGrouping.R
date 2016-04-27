#' Signal grouping
#'
#' Euclidean distances between m/z signals are hierarchically clustering using 
#' the average method and the composite spectrum groups determined by a absolute
#' error cutoff
#' 
#' @param spectrum.df a dataframe or matrix with two or more columns:
#' 1. Mass/ Mass-to-charge ratio
#' 2. Intensity
#' @param mzError interpeak absolute m/z error for signal grouping 
#' (Default = 0.001)
#' @return dataframe of m/z grouped signals, the m/z values of the input 
#' dataframe/ matrix peak groups are averaged and the signal intensities summed.
#' @export
signalGrouping <- function(spectrum.df=NULL, mzError=0.001, 
                           minPeaks = 5){
  # error handling
  if(is.null(spectrum.df)){
    stop("No spectrum matrix/dataframe supplied")    
  } else {
    hr <- fastcluster::hclust(dist(spectrum.df[, 1]), method = "median", members=NULL)
    # cut tree according to absolute m/z error
    spectrum_group <- cutree(hr, h=mzError)
     # calculate weighted mean of the m/z and sum signal within each peak group
    mass <- do.call(c, as.list(by(spectrum.df, as.factor(spectrum_group), function(x){
                                                   weighted.mean(x[, 1], x[, 2])})))
    grouped.df <- data.frame(mass = mass, 
                            intensity = tapply(spectrum.df[, 2], 
                                                 as.factor(spectrum_group), sum),
                              stringsAsFactors = F)
    #average any additional columns i.e. retention time
    if(ncol(spectrum.df) > 2){
    groupedCols <- apply(spectrum.df[, 3:ncol(spectrum.df), drop=F], 2, function(x) 
                         tapply(x, as.factor(spectrum_group), mean))
    grouped.df <- cbind(grouped.df, groupedCols)
    }
    if(nrow(grouped.df) > minPeaks){
      if(!is.null(colnames(spectrum.df))){
      colnames(grouped.df) <- colnames(spectrum.df)
      }
      return(grouped.df)
    } else {
      return("Less than minPeak")  
    }
  }
}
