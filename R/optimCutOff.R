#' optimum cutoff for a correlation or spectral similarity matrix
#' @param x matrix of correlation coefficients or spectral similarity values
#' column and row names must match
#' @param cutOffSeq numeric a vector of cut-off values to test. (default=seq(0.01, 1, 0.01) a numeric vector of length 100)
#' @param diffConsecVals numeric the scaled difference between consecutive values
#' to identify the plateau in the network density (default=1.5*10^-3 or 0.0015).
#' @return a list containing two named elements "estCutOff" a numeric estimated optimal cut-off value and "testData" a matrix
#' of the test results at each cut-off value. The function also plots the result.
#' @references Koh Aoki, Yoshiyuki Ogata, and Daisuke Shibata
#' Approaches for Extracting Practical Information from Gene Co-expression Networks in Plant Biology
#' Plant Cell Physiol (2007) 48 (3): 381-390 first published online January 23, 2007 doi:10.1093/pcp/pcm013
#' @export
optimCutOff <- function(x=NULL, cutOffSeq=seq(0.01, 1, 0.01), diffConsecVals=1.5*10^-3,
                        maxCutOff=0.95){
  if(!require(igraph)){
    stop('igraph package must be installed to use this function')
  }
  if(is.null(row.names(x)) | is.null(colnames(x))){
    stop('The matrix must have row and column names')
  }
  if(!all.equal(row.names(x), colnames(x))){
    stop('The row and column names must be equal')
  }
  # remove upper tri
  x[upper.tri(x, diag=TRUE)] <- 0
  # ID features above below corrThresh
  netData <- matrix(0, ncol=4, nrow=length(cutOffSeq))
  colnames(netData) <- c('nNodes', 'nEdges', 'nClust', 'netDense')
  row.names(netData) <- cutOffSeq
  for(i in 1:length(cutOffSeq)){
    arrIdxTmp <- which(abs(x) >= cutOffSeq[i], arr.ind = TRUE)
    arrIdxTmp <- cbind(colnames(x)[arrIdxTmp[, 1]], colnames(x)[arrIdxTmp[, 2]])
    if(nrow(arrIdxTmp) > 0){
      netTmp <- igraph::graph(as.vector(t(arrIdxTmp[, c(1, 2)])))
      netData[i, 'nNodes'] <- length(igraph::V(netTmp))
      netData[i, 'nEdges'] <- length(igraph::E(netTmp))
      netData[i, 'nClust'] <- igraph::clusters(netTmp)$no
    } 
  }
  # network density
  pConn <- {netData[, 'nNodes'] * {netData[, 'nNodes'] - 1}}/2
  netData[, 'netDense'] <- netData[, 'nEdges']/pConn
  netData[is.na(netData)] <- 0
  # estimate best cutoff value
  # identify plateau in the network density
  estCutOff <- 0
  multVal <- 1
  while(estCutOff == 0){
  diffConsecSeq <- abs(diff(netData[, 'netDense']/max(netData[, 'netDense'])))
  diffConsecIdx <- which(diffConsecSeq < {diffConsecVals * multVal} & diffConsecSeq > {diffConsecVals * {diffConsecVals * 0.1}} & as.numeric(names(diffConsecSeq)) > 0.6)
  firstDiffConsec <- ifelse(length(diffConsecIdx) == 0, NA, min(diffConsecIdx))
  if(!is.na(firstDiffConsec)){
  # at least two consecutive values to establish plateau
  firstDiffConsec <- ifelse({firstDiffConsec + 1} %in% diffConsecIdx, firstDiffConsec, NA)
  }
  multVal <- multVal + 0.01
  estCutOff <- ifelse(is.na(firstDiffConsec), 0, cutOffSeq[firstDiffConsec + 1])
    if(multVal > 2){
      estCutOff <- NA
      break
    }
  }
  if(is.na(estCutOff)){
    warning('Failed to find a plateau in the network density. Returning an approximation close to the maximum number of clusters on the upward slope.\n', immediate. = TRUE)
    flush.console()
    nClustTmp <- netData[, 'nClust']
    maxClustTmp <- max(nClustTmp)
    minClustTmp <- min(nClustTmp[cutOffSeq < 0.4])
    idxTmp <- c(max(which(nClustTmp == minClustTmp)), which(nClustTmp == maxClustTmp))
    medUpSlope <- round(quantile(idxTmp[1]:idxTmp[2], probs = seq(0, 1, 0.33))['66%'], 0)
    estCutOff <- cutOffSeq[medUpSlope]
  }
  if(estCutOff >= maxCutOff){
    warning('Estimated cut-off value greater than ', round(maxCutOff, 2), '. This indicates that a large proportion of the network nodes are highly related. In the case of spectral similarity for example this could indicate large numbers of similar spectra which were not removed during the combineMS2.removeContam step for example or biologically related such as lipids.\n', immediate. = TRUE)
    estCutOff <- maxCutOff
    }
  # estimate cutoff
  par(mfrow=c(2, 2))
  for(j in 1:ncol(netData)){
    plot(x=cutOffSeq, y=netData[, j], xlab='cut off value', pch=19, col='red', main=colnames(netData)[j], ylab=colnames(netData)[j])
    abline(v=rep(estCutOff, length(cutOffSeq)), col='blue')
    text(x=estCutOff, y=max(netData[, j])/2, paste0('estCutOff-', round(estCutOff, 2)), col='blue')  
  }
  par(mfrow=c(1, 1))
  return(list(estCutOff=estCutOff, testData=netData))
} # end function
