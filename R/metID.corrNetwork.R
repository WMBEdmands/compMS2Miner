#' CompMS2 correlation network generation
#' 
#' @description Uses the MS1features matched to MS2 data to generate a correlation network to view in compMS2explorer
#' 
#' @param object A "CompMS2" class object.  
#' @param peakTable a data.frame in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. 
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param corrThresh correlation coefficient threshold to group features within
#' a retention time cluster.
#' @param corrMethod character correlation method see \code{\link{cor}} for details. default "spearman".
#' @param delta numeric maximum p-value (following multiple testing correction) above 
#' which the null hypothesis (no correlation) is rejected.
#' @param MTC character Multiple Testing Correction default is "BH", see \code{\link{p.adjust.methods}} for
#' details of options. ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' Any p-values after multiple testing correction above the value of delta will have their
#' corresponding correlation coefficents replaced with zero.
#' @param maxNodes numeric above the maximum nodes the function will use the large graphing algorithm of \code{\link{igraph}}. See \code{\link{igraph::with_lgl}} else
#' the function uses the  Fruchterman-Reingold layout algorithm. See \code{\link{igraph::with_fr}} 
#' 
#' @return "CompMS2" class object with an additional network graph of any peakTable features above the correlation threshold.
#' 
#' @export
setGeneric("metID.corrNetwork", function(object, ...) standardGeneric("metID.corrNetwork"))

setMethod("metID.corrNetwork", signature = "CompMS2", function(object, peakTable=NULL, obsNames=NULL, 
                          corrThresh=0.6, corrMethod="spearman", delta=0.05, 
                          MTC="BH", maxNodes=300){
    # error handling
    stopifnot(!is.null(object))
    if(class(object) != "CompMS2"){
      stop('argument object must be a "CompMS2" class object')
    }
    if(is.null(obsNames)){
      stop('argument obsNames is missing with no default')
    } else if(is.null(peakTable)){
      stop('argument peakTable is missing. This should be the same peak-picker output table you initiated the CompMS2miner process with.')
    } 
    if(!require(igraph)){
      stop('the igraph package must be installed to calculate the correlation network. Try install.packages("igraph")')
    }
    # add parameters to object
    object@Parameters$corrThresh <- corrThresh
    object@Parameters$corrMethod <- corrMethod
    object@Parameters$deltaCorr <- delta
    object@Parameters$MTC <- MTC
    
    # match obsNames to peak table colnames
    obsIndx <- match(obsNames, colnames(peakTable))
    # if less than all matched then stop
    if(length(obsIndx) < length(obsNames)){
      stop(length(obsIndx), " of ", length(obsNames), 
           " observation names were matched in the peakTable column names, check the obsNames and peakTable column names")
    }
    # subset table
    obsTable <- peakTable[, obsIndx]
    
    if(any(is.na(obsTable))){
      message(sum(is.na(obsTable)), 'NAs in peakTable, replacing with zero prior to correlation matrix calculation.\n')
      flush.console()
      obsTable[is.na(obsTable)] <- 0
    }
    # calc correlation matrix
    message("Calculating correlation matrix for ", nrow(obsTable), " features\n\n")
    flush.console()
    cor.prob <- function(X, dfr = nrow(X) - 2, Method=corrMethod) {
      R <- cor(X, method=Method)
      r2 <- R^2
      Fstat <- r2 * dfr / (1 - r2)
      P <- 1 - pf(Fstat, 1, dfr)
      list(R, P)
    }
    # returns list R and P values
    corProbRes <- cor.prob(t(obsTable), Method=corrMethod)
    cor.m <- corProbRes[[1]]
    cor.m[apply(corProbRes[[2]], 2, p.adjust, 
                n=sum(lower.tri(cor.m)), method=MTC) > delta] <- 0
    colnames(cor.m)  <- peakTable[, 1]
    row.names(cor.m) <- peakTable[, 1]
    # replace upper tri with zero
    cor.m[upper.tri(cor.m, diag=T)] <- 0
    # ID features above below corrThresh
    sif <- apply(cor.m, 2, function(x){
      pos.tmp <- which(x > corrThresh)
      neg.tmp <- which(x < -corrThresh)
      if(length(pos.tmp) > 0){
        names(pos.tmp) <- x[pos.tmp]
      }
      if(length(neg.tmp) > 0){
        names(neg.tmp) <- x[neg.tmp]
      }
      return(list(pos=pos.tmp, neg=neg.tmp))
    })
    # melt list result
    sif.df <- reshape2::melt(sif)
    # add correlation value
    sif.df[, 4] <- as.numeric(gsub('.+pos\\.|.+neg\\.', '', names(unlist(sif))))
    # create sif file names
    sif.df[, 5] <- names(sif)[sif.df[, 1]]
    
    # add in features with no correlation above thresh
    # noCorrFeat <- setdiff(peakTable[, 1], unique(c(sif.df[, 3], sif.df[, 5])))
    # noCorrFeatM <- matrix(0, nrow=length(noCorrFeat), ncol=5)
    # colnames(noCorrFeatM) <- colnames(sif.df)
    # noCorrFeatM[, 1] <- noCorrFeat
    # noCorrFeatM[, 3] <- noCorrFeat
    # sif.df <- rbind(sif.df, noCorrFeatM)
    
    netTmp <- igraph::graph(as.vector(t(sif.df[, c(1, 3)])))
    nNodes <- length(igraph::V(netTmp))
    if(nNodes <= maxNodes){
      message('less than ', maxNodes, ' nodes using Fruchterman-Reingold layout. see ?igraph::with_fr()\n')
      flush.console()
      layoutTmp <- igraph::layout_(netTmp, with_fr()) 
    } else {
      message('more than ', maxNodes, ' nodes using large-layout function. see ?igraph::with_lgl()\n')
      flush.console()
      layoutTmp <- igraph::layout_(netTmp, with_lgl())   
    }
    message(nNodes, 
            " nodes with ", length(igraph::E(netTmp)), " edges identified at a correlation threshold >= ", corrThresh, " (", corrMethod, ', p/q value <= ', delta, ', MTC: ', MTC, ')')
    flush.console()
    
    object@network <- list(networkGraph=netTmp, layout=layoutTmp)
    return(object)
   
}) # end function