#' remove any possible contaminants identified as repeating sequences of isobaric ions of 
#' high spectral similarity.
#' 
#' @details The function produces plots to visualize the contaminants identified.
#' If you suspect all isobaric ions across the gradient are contaminants then
#' you can set the argument maxRtGap to Infinite (maxRtGap=Inf). N.B. Any spectrum which is not sufficiently similar (specSimScore) to the rest
#' of the potential contaminants across the gradient will not be erroneously removed.
#' @param object a "compMS2" class object.
#' @param ms1Abs numeric ms1 mass to charge absolute error for hierarchical 
#' clustering of MS1 masses (default = 0.01). Utilizes the median method of the
#' \code{\link{hclust}} function of the fastcluster package.
#' @param maxRtGap numeric maximum retention time gap (in seconds) between two 
#' isobaric ions. Sequences of possible contaminants will be identified using this
#' difference in retention time (default = 60). 
#' @param nContams numeric number of isobaric ions in a sequence separated by a maximum
#' retention time gap (maxRtGap). Any sequences of isobaric ions greater than or equal 
#' to this number will be removed (default = 10).
#' @param minSimScore numeric minimum spectral similarity score (values between 0-1).
#' If any isobar in a sequence of possible contaminants is below this minimum 
#' mean dot product similarity score
#' then it will not be removed. This is to ensure that only true isobaric contaminants
#' are removed and spectra which have been grouped amongst them are not erroneously
#' removed (default = 0.8).
#' @param remContam logical should possible contaminant spectra be automatically
#' removed from the object (default = TRUE), If FALSE the contaminant plot
#' will still be printed but the spectra will not be removed. In this way the user
#' can interactively determine suitable parameters prior to spectrum removal.
#' @examples 
#' compMS2contamRem <- combineMS2(compMS2Example, 'removeContam', maxRtGap=Inf, 
#'                                nContams=4)
#' @return a "compMS2" class object with the composite spectra and all metID 
#' information of any contaminants identified removed. Any correlation or spectral
#' similarity networks will also have to be recalculated.
#' @export
setGeneric("combineMS2.removeContam", function(object, ...) standardGeneric("combineMS2.removeContam"))

setMethod("combineMS2.removeContam", signature = "compMS2", 
          function(object, ms1Abs=0.01, maxRtGap=60, nContams=10,
                   minSimScore=0.8, remContam=TRUE){
  # error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  message('Identifying possible contaminants (sequences of isobaric ions with high spectral similarity) using the following parameters:\n', '1. absolute MS1 mass error: ', 
          round(ms1Abs, 3), '\n2. maximum retention time gap between isobar: ', 
          round(maxRtGap, 0), ' (seconds)', 
          '\n3. minimum mean similarity score (dot product): ', 
          round(minSimScore, 1), '\n\n')
  flush.console()
rtsTmp <- sapply(metaData(object), function(x) x[grep('_MS1_RT', names(x))][[1]][1])
mzsTmp <- sapply(metaData(object), function(x) x[grep('_MS1_mz', names(x))][[1]][1])
mzsTmp <- mzsTmp[order(rtsTmp)]
rtsTmp <- rtsTmp[order(rtsTmp)]

plot(rtsTmp, mzsTmp, ylab='m/z', xlab='retentionTime', pch=19)
hr <- fastcluster::hclust(dist(mzsTmp), method = "median", members=NULL)

# cut tree according to absolute m/z error
massGroup <- cutree(hr, h=ms1Abs)
freqGroup <- table(massGroup)
possContams <- tapply(rtsTmp, massGroup, function(x){
 spacingTmp <- rle(diff(x) < maxRtGap)
 possContamsTmp <- spacingTmp$lengths >= (nContams - 1) & spacingTmp$values
 if(any(possContamsTmp)){
     cumSumSpacL <- cumsum(spacingTmp$lengths) + 1
     stopTmp <- cumSumSpacL[possContamsTmp]
     startTmp <- c(1, cumSumSpacL)[which(possContamsTmp)]
     contamsTmp <- do.call(c, mapply(seq, from=startTmp, to=stopTmp, by=1, SIMPLIFY = FALSE))
     specNamesTmp <- names(x)[contamsTmp]
     allSpecTmp <- do.call(rbind, compSpectra(object)[specNamesTmp])
     allSpecTmp <- as.matrix(allSpecTmp[, c('mass', 'intensity'), drop=FALSE])
     specNamesVecTmp <- gsub("\\.[^\\.]*$", '', row.names(allSpecTmp))
     # check spectral similarity
     dpScores <- suppressMessages(dotProdMatrix(allSpecTmp, specNamesVecTmp, binSizeMS2 = 0.1))
     # column means
     dpScores <- colMeans(dpScores)
     specNamesTmp <- specNamesTmp[dpScores >= minSimScore]
     if(length(specNamesTmp) > 0){
     return(specNamesTmp)
     }
 }})
possContams <- unlist(possContams)
if(length(possContams) > 0){
  message('Of ', length(mzsTmp), ' composite spectra.\n', length(possContams),
          ' possible contaminants (', round(length(possContams)/length(mzsTmp), 1) * 100, '%) were detected.\n', 
          '(i.e. greater than or equal to ', nContams, ' isobaric ions in a sequence separated by no greater than ', maxRtGap, ' seconds)\n\n')
  flush.console()
  
  possContamIndx <- names(massGroup) %in% possContams
    
  points(rtsTmp, mzsTmp, pch=19, col=ifelse(possContamIndx, 'red', 'black'))
  legend('topleft', c('possible contaminant'), pch=19, col='red')
  
  # message('CAUTION: do you wish to permanently remove the possible contaminants',
  #         ' from the "compMS2" class object and therefore from further downstream analysis?\n', ifelse(length(network(object)), 'You may also have to recalculate the correlation and spectral similarity networks.\n', ''), 'This cannot be undone.\ntype (Y/N) and press [enter] to continue:')
  # flush.console()
  # remContams <- readline()
  # if(remContams == 'Y'){
  if(remContam == TRUE){
    message('Removing possible contaminant spectra from the "compMS2" class object and therefore from further downstream analysis\n')
    flush.console()
  # index of contaminant spectra
  specKeep <- setdiff(names(compSpectra(object)), possContams)
  # selectively keep only non-contaminant elements
  object <- subsetCompMS2(object, specKeep)
  }
  # add details of contaminants to comments of any remaining mass groups
  massGroup <- massGroup[possContamIndx == FALSE]
  freqGroup <- table(massGroup)
  freqGroup <- freqGroup[freqGroup >= nContams]
  if(length(freqGroup) > 0){
  massGroup <- massGroup[massGroup %in% names(freqGroup)]  
  # assemble labels
  possContamLabel <- sapply(massGroup, function(x){ 
    redunSpec <- setdiff(names(massGroup)[massGroup == x], names(x))
    paste0(length(redunSpec), ' possible contaminants. Isomer of: ', 
           paste0(redunSpec, collapse = ', '))})
  metIDcomments <- Comments(object)
  indxTmp <- match(metIDcomments$compSpectrum, names(massGroup))
  metIDcomments$user_comments[!is.na(indxTmp)] <- possContamLabel[indxTmp[!is.na(indxTmp)]]
  Comments(object) <- metIDcomments
  }
} else {
  message('No potential contaminants were detected.\n')
  flush.console()
}

return(object)
}) # end function
