#' CompMS2 correlation network generation
#' 
#' @description Uses the MS1features matched to MS2 data to generate a correlation network to view in compMS2Explorer
#' 
#' @param object A "compMS2" class object.  
#' @param peakTable a data.frame in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. The first 3 columns must consist of:
#' \enumerate{
#'  \item EIC number or unique peak identifier.
#'  \item mass-to-charge ratio of peak group.
#'  \item median/ peak apex retention time in seconds. 
#'  }
#' These columns are utilized in the final network visualization.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param corrThresh correlation coefficient threshold to group features within
#' a retention time cluster. If no value is supplied the default is to estimate
#' an optimal cut-off value based on network attributes from a series of correlation cut-off values. This is based upon the relationship between the 
#' number of nodes, edges, number of clusters and the network density (i.e. n actual edges/n potential edges) at each correlation cut-off value. A plot showing
#' the result of this estimation will be generated.
#' @param corrMethod character correlation method see \code{\link{cor}} for details. default "spearman".
#' @param delta numeric maximum p-value (following multiple testing correction) above #' which the null hypothesis (no correlation) is rejected.
#' @param MTC character Multiple Testing Correction default is "none", see \code{\link{p.adjust.methods}} for
#' details of options. ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' Any p-values after multiple testing correction above the value of delta will have their
#' corresponding correlation coefficents replaced with zero.
#' @param maxNodes numeric above the maximum nodes the function will use the large graphing algorithm of \code{\link{igraph}}. See \code{\link{igraph::with_lgl}} else
#' the function uses the  Fruchterman-Reingold layout algorithm. See \code{\link{igraph::with_fr}} 
#' @param MS2only numeric 3 options (1-3) if 1 All EICs above corrThresh returned, if 2 only non-MS2 matched EICs which are first neighbours of at least one MS2-matched MS2 are returned, if 3 only MS2-matched EICs are returned (default = 1).
#' @param minClustSize numeric minimum number of connected nodes that is cluster size (default =3). If a cluster
#' of nodes is less than this number then it will be removed.
#' @return "compMS2" class object with an additional network graph of any peakTable features above the correlation threshold.
#' @references Koh Aoki, Yoshiyuki Ogata, and Daisuke Shibata
#' Approaches for Extracting Practical Information from Gene Co-expression Networks in Plant Biology
#' Plant Cell Physiol (2007) 48 (3): 381-390 first published online January 23, 2007 doi:10.1093/pcp/pcm013
#' @export
setGeneric("metID.corrNetwork", function(object, ...) standardGeneric("metID.corrNetwork"))

setMethod("metID.corrNetwork", signature = "compMS2", function(object, peakTable=NULL, obsNames=NULL, corrThresh=NULL, corrMethod="spearman", delta=0.05, MTC="none",
maxNodes=300, MS2only=1, minClustSize=2){
    # error handling
    stopifnot(!is.null(object))
    if(class(object) != "compMS2"){
      stop('argument object must be a "compMS2" class object')
    }
    if(is.null(obsNames)){
      stop('argument obsNames is missing with no default')
    } else if(is.null(peakTable)){
      stop('argument peakTable is missing. This should be the same peak-picker output table you initiated the compMS2Miner process with.')
    } 
    if(!require(igraph)){
      stop('the igraph package must be installed to calculate the correlation network. Try install.packages("igraph")')
    }
    if((MS2only %in% c(1:3)) == FALSE){
      stop('The MS2only argument must be one of 3 options: 1, 2 or 3.')
    }
    # add parameters to object
    Parameters(object)$corrMethod <- corrMethod
    Parameters(object)$deltaCorr <- delta
    Parameters(object)$MTC <- MTC
    # index ms2 matched
    ms2MatchedEICs <- as.numeric(gsub('.+_', '', names(compSpectra(object))))
    
    if(MS2only == 3){
      # subset to include MS2 matched
      # subset peakTable to only include MS2 matched
      peakTable <- peakTable[peakTable$EICno %in% ms2MatchedEICs, ]
    }
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
    message("Calculating correlation matrix for ", 
            prettyNum(nrow(obsTable), big.mark = ','), " features\n\n")
    flush.console()
    if(nrow(obsTable) > 10000){
      warning('This is a large matrix and the correlation matrix may take some time or your R session may run out of memory.\n', immediate.=TRUE)
      flush.console()
    }
    # from https://stat.ethz.ch/pipermail/r-help/2000-January/009758.html
    cor.prob <- function(X, dfr = nrow(X) - 2, Method=corrMethod) {
      R <- cor(X, method=Method)
      r2 <- R^2
      Fstat <- r2 * dfr / (1 - r2)
      P <- 1 - pf(Fstat, 1, dfr)
      list(R, P)
    }
    # END from https://stat.ethz.ch/pipermail/r-help/2000-January/009758.html
    # returns list R and P values
    corProbRes <- cor.prob(t(obsTable), Method=corrMethod)
    cor.m <- corProbRes[[1]]
    cor.m[apply(corProbRes[[2]], 2, p.adjust, 
                n=sum(lower.tri(cor.m)), method=MTC) > delta] <- 0
    colnames(cor.m)  <- peakTable[, 1]
    row.names(cor.m) <- peakTable[, 1]
    cor.m[is.na(cor.m)] <- 0
    # identify optimum correlation cutoff
    if(is.null(corrThresh)){
      cat('Estimating optimal correlation coefficient based on network attributes from a sequence of cut-off values.\n')
      if(ncol(cor.m) > 2000){
        cat('\nLarge correlation matrix', prettyNum(ncol(cor.m), big.mark = ','), 
            'x', prettyNum(ncol(cor.m), big.mark = ','), 
            ': starting from a cutoff of 0.5 correlation coefficient.\n')
        cutOffSeqTmp <- seq(0.5, 1, 0.05)
      } else {
        cutOffSeqTmp <- seq(0.01, 1, 0.01)
      }
      corrThresh <- optimCutOff(cor.m, cutOffSeq = cutOffSeqTmp)$estCutOff
      cat('\nestimate of optimal corrThresh value:', corrThresh, '.\n')
    }
    Parameters(object)$corrThresh <- corrThresh
    # replace upper tri with zero
    cor.m[upper.tri(cor.m, diag=TRUE)] <- 0
    # ID features above below corrThresh
    posCorrTmp <- which(cor.m > corrThresh, arr.ind = TRUE)
    negCorrTmp <- which(cor.m < {-corrThresh}, arr.ind = TRUE)
    if(nrow(posCorrTmp) == 0 & nrow(negCorrTmp) == 0){
      stop('no correlations above a threshold of ', corrThresh)
    }
    # correlation coefficient value
    posCorrTmp <- cbind(posCorrTmp, cor.m[posCorrTmp])
    posCorrTmp <- cbind(posCorrTmp, 'pos')
    if(nrow(negCorrTmp) > 0){
    # correlation coefficient value
    negCorrTmp <- cbind(negCorrTmp, cor.m[negCorrTmp])
    negCorrTmp <- cbind(negCorrTmp, 'neg')
    sif.df <- rbind(posCorrTmp, negCorrTmp)
    } else {
    sif.df <- posCorrTmp
    }
    
    sif.df[, 1] <- peakTable[as.numeric(sif.df[, 1]), 1]
    sif.df[, 2] <- peakTable[as.numeric(sif.df[, 2]), 1]
    # sif <- apply(cor.m, 2, function(x){
    #   pos.tmp <- which(x > corrThresh)
    #   neg.tmp <- which(x < -corrThresh)
    #   if(length(pos.tmp) > 0){
    #     names(pos.tmp) <- x[pos.tmp]
    #   }
    #   if(length(neg.tmp) > 0){
    #     names(neg.tmp) <- x[neg.tmp]
    #   }
    #   return(list(pos=pos.tmp, neg=neg.tmp))
    # })
    # # melt list result
    # sif.df <- reshape2::melt(sif)
    # add correlation value
    # create sif file names
    # sif.df[, 1] <- colnames(cor.m)[as.numeric(sif.df[, 1]]
    
    # add average intensity
    avInt <- unique(as.vector(sif.df[, 1:2]))
    avInt <- rowMeans(peakTable[peakTable[, 1] %in% avInt, obsNames])
    names(avInt) <- ifelse(names(avInt) %in% ms2MatchedEICs, 
                           paste0('CC_', names(avInt)), 
                           paste0('noMS2_', names(avInt)))
    # add prefix to names
    indxTmp <- sif.df[, 1] %in% ms2MatchedEICs
    sif.df[, 1] <- ifelse(indxTmp, paste0('CC_', sif.df[, 1]), 
                          paste0('noMS2_', sif.df[, 1]))
    indxTmp <- sif.df[, 2] %in% ms2MatchedEICs
    sif.df[, 2] <- ifelse(indxTmp, paste0('CC_', sif.df[, 2]), 
                          paste0('noMS2_', sif.df[, 2]))
    
    # add in features with no correlation above thresh
    # noCorrFeat <- setdiff(peakTable[, 1], unique(c(sif.df[, 3], sif.df[, 5])))
    # noCorrFeatM <- matrix(0, nrow=length(noCorrFeat), ncol=5)
    # colnames(noCorrFeatM) <- colnames(sif.df)
    # noCorrFeatM[, 1] <- noCorrFeat
    # noCorrFeatM[, 3] <- noCorrFeat
    # sif.df <- rbind(sif.df, noCorrFeatM)
    
    netTmp <- igraph::graph(as.vector(t(sif.df[, c(1, 2)])))
    # add corr coeff
    igraph::E(netTmp)$value <- sif.df[, 3]
    # add average intensity
    idxTmp <- match(names(igraph::V(netTmp)), names(avInt))
    igraph::V(netTmp)$avInt <- avInt[idxTmp]
    # neighbours of MS2 matched
    if(MS2only == 2){
    eicsMS2 <- which(igraph::V(netTmp)$name %in% names(compSpectra(object)))
    adjVertTmp <- c(eicsMS2, unlist(adjacent_vertices(netTmp, eicsMS2, mode='all')))
    netTmp <- igraph::induced_subgraph(graph=netTmp, vids=unique(adjVertTmp))
    }
    
    # if necessary remove clusters
    if(minClustSize > 2){
      clustsTmp <- igraph::clusters(netTmp)
      beforeClust <- length(clustsTmp$csize)
      indxTmp <- which(clustsTmp$csize >= minClustSize)
      retVert <- names(clustsTmp$membership)[clustsTmp$membership %in% indxTmp]
      # take subgraph
      netTmp <- igraph::induced_subgraph(graph=netTmp, vids=retVert)
      afterRem <- length(igraph::clusters(netTmp)$csize)
      message(beforeClust - afterRem, ' clusters out of ', beforeClust, ' consisted of less than ', minClustSize, ' nodes and were removed.\n')
      flush.console()
    }
    
    # add any noMS2 matched to object
    if(MS2only %in% c(1, 2)){
      # add any noMS2 matched to object
      message('Adding the non-MS2 matched EICs to the "compMS2" object.\n')
      flush.console()
      if(!is.null(network(object)$specSimGraph)){
        warning('The spectral similarity network will have to be re-calculated.\n',
                immediate. = TRUE)
        flush.console()
      }
      # 1. subset to remove any noMS2 already present (i.e. a previous corrNetwork or specSim network)
      namesSub <- names(compSpectra(object))[grep('^CC_', names(compSpectra(object)))]
      object <- subsetCompMS2(object, namesSub)
      # 2. add noMS2 names to metaData, compSpectra, DBanno
      namesTmp <- igraph::V(netTmp)$name
      namesTmp <- namesTmp[grep('^noMS2', namesTmp)]
      eicNos <- gsub('^noMS2_', '', namesTmp)
      eicMzRt <- peakTable[match(eicNos, peakTable[, 1]), 1:4, drop=FALSE]
      # see if the 4 th column contains adduct information
      adducts <- any(grepl('M\\+|M\\-', peakTable[, 4]))
      if(adducts == FALSE){
      eicMzRt[, 4] <- ''
      }
      object <- addNoMS2(object, namesTmp, eicMzRt)
    }
   
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
    # add EIC no., mass and RT to layout
    indxTmp <- match(as.numeric(gsub('.+_', '', igraph::V(netTmp)$name)), peakTable[, 1])
    layoutTmp <- cbind(layoutTmp, as.matrix(peakTable[indxTmp, 1:3]))
    network(object)$corrNetworkGraph <- netTmp
    network(object)$corrLayout <- data.frame(layoutTmp, stringsAsFactors = FALSE)
    
    return(object)
}) # end function
