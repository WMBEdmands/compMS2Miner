#' dot product matrix calculation
#' @param allSpectra a numeric matrix consisting of two columns 1. mass and 2. intensity
#' @param spectraNames character names of individual spectra to compare must equal number of rows of allSpectra
#' @param binSizeMS2 numeric the MS2 bin size to bin MS2 data prior to dot product calculation (default = 0.1 Da).
#' @return a matrix of equal dimension corresponding to the number of factor levels f
dotProdMatrix <- function(allSpectra=NULL, spectraNames=NULL, binSizeMS2=NULL){
  # error handling
  stopifnot(is.matrix(allSpectra))
  stopifnot(is.character(spectraNames))
  message('Calculating dot product matrix ', length(unique(spectraNames)), ' spectra\n')
  flush.console()
  # padded integer labels
  labelsTmp <- paste0('(', seq(binSizeMS2, (1000 - binSizeMS2), binSizeMS2), ',', seq((2 * binSizeMS2), 1000, binSizeMS2), ']')
  massBinsIndivTmp <- cut(allSpectra[, 1], breaks=seq(binSizeMS2, 1000, binSizeMS2), labels=labelsTmp)   
  # empty bins
  indivSpecVec <- tapply(allSpectra[, 2], paste0(spectraNames, massBinsIndivTmp), sum)
  # identify any absent bins
  allBinNames <- paste0(rep(unique(spectraNames), each=length(labelsTmp)), rep(labelsTmp, length(unique(spectraNames))))
  # add absent bins as zeros
  allBinsTmp <- rep(0, length(allBinNames))
  names(allBinsTmp) <- allBinNames
  # ensure indivSpecVec is in right order
  allBinsTmp[match(names(indivSpecVec), allBinNames)] <- indivSpecVec
  
  indivSpecMat <- matrix(allBinsTmp, byrow=F, nrow=length(labelsTmp))
  # mean all pairwise dotproducts
  dotProdMat <- crossprod(indivSpecMat)
  sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), nrow=nrow(dotProdMat), 
                          ncol=ncol(dotProdMat), byrow = T) 
  
  dotProdsTmp <- dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
  row.names(dotProdsTmp) <- unique(spectraNames)
  colnames(dotProdsTmp) <- unique(spectraNames)
  return(dotProdsTmp)
} # end function