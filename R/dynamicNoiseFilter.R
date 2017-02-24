#' Dynamic Noise filtration 
#' 
#' @param spectrum.df a dataframe or matrix with two columns:
#' 1. Mass/ Mass-to-charge ratio
#' 2. Intensity
#' @param DNF dynamic noise filter minimum signal to noise threshold 
#' (default = 2), calculated as the ratio between the linear model predicted 
#' intensity value and the actual intensity.
#' @param minPeaks minimum number of signal peaks following dynamic 
#' noise filtration (default = 5).
#' @param maxPeaks maximum number of signal peaks the function will continue
#' until both the minimum DNF signal to noise ratio is exceeding and the number
#' of peaks is lower than the maximum (default = 5).
#' 
#' @return a list containing 3 objects:
#' \enumerate{
#' \item Above.noise The dynamic noise filtered matrix/ dataframe 
#' \item metaData a dataframe with the following column names:
#'        1. Noise.level the noise level determined by the dynamic noise filter 
#'           function.
#'        2. IntCompSpec Total intensity composite spectrum.
#'        3. TotalIntSNR Sparse ion signal to noise ratio 
#'        (mean intensity/ stdev intensity)
#'        4. nPeaks number of peaks in composite spectrum
#' \item aboveMinPeaks Logical are the number of signals above the minimum level}
#' @details  Dynamic noise filter adapted from the method described in Xu H. and 
#' Frietas M. "A Dynamic Noise Level Algorithm for Spectral Screening of 
#' Peptide MS/MS Spectra" 2010 BMC Bioinformatics. 
#' 
#' The function iteratively calculates linear models starting from 
#' the median value of the lower half of all intensities in the spectrum.df. 
#' The linear model is used to predict the next peak intensity and ratio is 
#' calculated between the predicted and actual intensity value. 
#' 
#' Assuming that all preceeding intensities included in the linear model 
#' are noise, the signal to noise ratio between the predicted and actual values 
#' should exceed the minimum signal to noise ratio (default DNF = 2). 
#' 
#' The function continues until either the DNF value minimum has been exceeded 
#' and is also below the maxPeaks or maximum number of peaks value. As the 
#' function must necessarily calculate potentially hundreds of linear models the 
#' RcppEigen package is used to increase the speed of computation.
#' 
#' @export
dynamicNoiseFilter <- function(spectrum.df=NULL, DNF=2, minPeaks=5, 
                               maxPeaks=20, minInt=100){
  # error handling
  if(is.null(spectrum.df)){
    stop("No spectrum matrix/dataframe supplied")    
  } else {
    # rank matrix/ dataframe by intensity
    intOrder<-order(spectrum.df[, 2])
    spectrum.df <- spectrum.df[intOrder, , drop=FALSE]
    # median bottom half of intensity values
#     medBottomHalf <- median(head(spectrum.df[, 2],
#                                  n=nrow(spectrum.df)/2))
#     medBottomHalf <- which(spectrum.df[, 2] >= medBottomHalf)[1]
    minIntIndx <- which(spectrum.df[, 2] >= minInt)[1]
    peakIndx <- seq(1, nrow(spectrum.df), 1)
    minIntIndx <- ifelse(is.na(minIntIndx), nrow(spectrum.df), minIntIndx)
    minIntIndx <- ifelse(minIntIndx == 1, 2, minIntIndx)
    # break loop if higher DNF and also less than maxPeaks
#     for(k in medBottomHalf:(nrow(spectrum.df)-1)){
if(minIntIndx < (nrow(spectrum.df)-1)){
for(k in minIntIndx:(nrow(spectrum.df)-1)){
      # calc linear model rcppeigen 
      fit <- coef(RcppEigen::fastLm(as.numeric(spectrum.df[1:k, 2]) 
                                    ~ peakIndx[1:k]))
      # predicted intensity model from intercept
      PredInt <- fit[1]+(fit[2])*(k+1)
      # calc Signal to noise ratio predicted vs. actual
      SNR <- spectrum.df[k+1, 2]/PredInt
      # if SNR reached and below max number of peaks break loop
      if(SNR >= DNF & nrow(spectrum.df) - (k+1) < maxPeaks){
        Noise.level <- as.numeric(spectrum.df[k+1, 2])
        break
      } else {
        Noise.level <- as.numeric(spectrum.df[k+1, 2])
      }
    }
} else {
  Noise.level <- spectrum.df[minIntIndx, 2]
}
    # filter by DNF noise filter level
    Noise.indx <- which(spectrum.df[, 2] >= Noise.level)
    spectrum.df <- spectrum.df[Noise.indx, , drop=FALSE]
    # sort by m/z
    spectrum.df <- spectrum.df[order(spectrum.df[,1]), , drop=FALSE]
    # number of peaks higher than minimum
    aboveMinPeaks <- nrow(spectrum.df) >= minPeaks
    # Total intensity composite spectrum 
    IntCompSpec <- sum(spectrum.df[, 2])
    # Sparse ion signal to noise ratio (mean intensity/ stdev intensity)
    if(nrow(spectrum.df) <= 1){
      TotalIntSNR <- 0
    } else {
      TotalIntSNR <- mean(spectrum.df[, 2], sd(spectrum.df[, 2]))
    }
    DNF.tmp<-list(Above.noise=spectrum.df, 
                  metaData=data.frame(Noise.level=Noise.level, 
                                      IntCompSpec=IntCompSpec,
                                      TotalIntSNR=TotalIntSNR,
                                      nPeaks=nrow(spectrum.df), 
                                      stringsAsFactors = FALSE),
                  aboveMinPeaks=aboveMinPeaks)
    return(DNF.tmp)
  }  
}
