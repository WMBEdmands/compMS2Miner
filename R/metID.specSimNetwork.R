#' CompMS2 spectral similarity network generation
#' 
#' @description generates a dot product spectral similarity network from the MS/MS fragmentation data. The resulting spectral similarity network can then be viewed in compMS2explorer. Utilizes the \code{\link{graph}} function of the \code{\link{igraph}} package.
#' 
#' @param object A "CompMS2" class object.  
#' @param minDotProdThresh minimum dot product spectral similarity score (default = 0.8).
#' @return "CompMS2" class object with an additional network graph of any peakTable features above the correlation threshold.
#' 
#' @export
setGeneric("metID.specSimNetwork", function(object, ...) standardGeneric("metID.specSimNetwork"))

setMethod("metID.specSimNetwork", signature = "CompMS2", function(object, minDotProdThresh=0.8){
    # error handling
    stopifnot(!is.null(object))
    if(class(object) != "CompMS2"){
      stop('argument object must be a "CompMS2" class object')
    }
  if(!require(igraph)){
    stop('the package igraph must be installed to utilize this function.')
  }
  # add parameters to object
  object@Parameters$minDotProdThresh <- minDotProdThresh
  
  # all constituent spectra
  allSpecTmp <- do.call(rbind, lapply(object@compSpectra, function(x){
   massIntTmp <- x[, c('mass', 'intensity')]
   #  preFragIntTmp <- x[, c('Precursorfrag.diff', 'intensity')]
   #  preFragIntTmp[, 1] <- as.numeric(preFragIntTmp[, 1])
   # preFragIntTmp <- preFragIntTmp[preFragIntTmp[, 1] > 2, , drop=F]
   # # colnames(preFragIntTmp)[1] <- 'mass'
   #  return(preFragIntTmp)
   }))
   specNamesVecTmp <- gsub('\\..+', '', row.names(allSpecTmp))
    # allSpecTmp <- cbind(allSpecTmp, specNamesVecTmp)
    # padded integer labels
    labelsTmp <- paste0(sprintf("(%04d", 1:999), ',', sprintf("%04d", 2:1000), ']')
    massBinsIndivTmp <- cut(allSpecTmp[, 1], breaks=seq(1, 1000, 1), labels=labelsTmp)
    # empty bins
    emptyBins <- cut(seq(2, 1000, 1), breaks=seq(1, 1000, 1), labels=labelsTmp)
    indivSpecVec <- tapply(allSpecTmp[, 2], paste0(specNamesVecTmp, massBinsIndivTmp), sum)
    # identify any absent bins
    allBinNames <- paste0(rep(unique(specNamesVecTmp), each=length(emptyBins)), rep(emptyBins, length(unique(specNamesVecTmp))))
    # add absent bins as zeros
    allBinsTmp <- rep(0, length(allBinNames))
    names(allBinsTmp) <- allBinNames
    allBinsTmp[which(allBinNames %in% names(indivSpecVec))] <- indivSpecVec
    indivSpecMat <- matrix(allBinsTmp, byrow=F, nrow=length(emptyBins))
    # mean all pairwise dotproducts
    # dotProdMat <- t(indivSpecMat) %*% indivSpecMat
    dotProdMat <- crossprod(indivSpecMat)
    sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), nrow=nrow(dotProdMat), 
                            ncol=ncol(dotProdMat), byrow = T) 
    
    fragsDotProds <- dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
    
    colnames(fragsDotProds)  <- names(object@compSpectra)
    row.names(fragsDotProds) <- names(object@compSpectra)
    
    # 2. precursor - neutral losses
    
    # all constituent spectra
    allSpecTmp <- do.call(rbind, lapply(object@compSpectra, function(x){
      #massIntTmp <- x[, c('mass', 'intensity')]
      preFragIntTmp <- x[, c('Precursorfrag.diff', 'intensity'), drop=F]
      preFragIntTmp[, 1] <- as.numeric(preFragIntTmp[, 1])
      preFragIntTmp <- preFragIntTmp[preFragIntTmp[, 1] > 2, , drop=F]
      if(nrow(preFragIntTmp) == 0){
      preFragIntTmp <- data.frame(Precursorfrag.diff=2, intensity=1, stringsAsFactors = F)  
      }
      # preFragIntTmp <- cbind(preFragIntTmp, name=names(object@compSpectra)[x])
      return(preFragIntTmp)
      # }
      # colnames(preFragIntTmp)[1] <- 'mass'
    }))
    specNamesVecTmp <- gsub('\\..+', '', row.names(allSpecTmp))
    # allSpecTmp <- cbind(allSpecTmp, specNamesVecTmp)
    # padded integer labels
    labelsTmp <- paste0(sprintf("(%04d", 1:999), ',', sprintf("%04d", 2:1000), ']')
    massBinsIndivTmp <- cut(allSpecTmp[, 1], breaks=seq(1, 1000, 1), labels=labelsTmp)
    # empty bins
    emptyBins <- cut(seq(2, 1000, 1), breaks=seq(1, 1000, 1), labels=labelsTmp)
    indivSpecVec <- tapply(allSpecTmp[, 2], paste0(specNamesVecTmp, massBinsIndivTmp), sum)
    # identify any absent bins
    allBinNames <- paste0(rep(unique(specNamesVecTmp), each=length(emptyBins)), rep(emptyBins, length(unique(specNamesVecTmp))))
    # add absent bins as zeros
    allBinsTmp <- rep(0, length(allBinNames))
    names(allBinsTmp) <- allBinNames
    allBinsTmp[which(allBinNames %in% names(indivSpecVec))] <- indivSpecVec
    indivSpecMat <- matrix(allBinsTmp, byrow=F, nrow=length(emptyBins))
    # mean all pairwise dotproducts
    # dotProdMat <- t(indivSpecMat) %*% indivSpecMat
    dotProdMat <- crossprod(indivSpecMat)
    sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), nrow=nrow(dotProdMat), 
                            ncol=ncol(dotProdMat), byrow = T) 
    
    preFragDotProds <- dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
    
    colnames(preFragDotProds)  <- names(object@compSpectra)
    row.names(preFragDotProds) <- names(object@compSpectra)
    
    # replace upper tri with zero
    fragsDotProds[upper.tri(fragsDotProds, diag=T)] <- 0
    # replace upper tri with zero
    preFragDotProds[upper.tri(preFragDotProds, diag=T)] <- 0
    
    # which dot product is larger fragment similarity or neutral loss?
    indxTmp <- preFragDotProds > fragsDotProds 
    # replace values, remove negatives and add prefragdotprods necessary and make 
    # prefrag dot products negative to distinguish them.
    fragsDotProds[fragsDotProds < 0] <- 0
    fragsDotProds[indxTmp] <- -preFragDotProds[indxTmp]
    
    # make preFragDot
    # calculate pca model
    # pcDp <- pcaMethods::pca(preFragDotProds, method="svd", nPcs=2, cv='q2', center = F)
    
    # # ID features above below corrThresh
    sif <- apply(fragsDotProds, 2, function(x){
      fragsTmp <- which(x >= minDotProdThresh)
      neutTmp <- which(x <= -minDotProdThresh)
      if(length(fragsTmp) > 0){
        names(fragsTmp) <- x[fragsTmp]
      }
      if(length(neutTmp) > 0){
        names(neutTmp) <- x[neutTmp]
      }
      return(list(frags=fragsTmp, neutLoss=neutTmp))
    })
    # melt list result
    sif.df <- reshape2::melt(sif)
    # add correlation value
    sif.df[, 4] <- as.numeric(gsub('.+frags\\.|.+neutLoss\\.', '', names(unlist(sif))))
    # create sif file names
    sif.df[, 5] <- names(sif)[sif.df[, 1]]
    netTmp <- igraph::graph(as.vector(t(sif.df[, c(3, 5)])))
    # add edge colours to igraph object based on fragment or neutral loss
    igraph::E(netTmp)$color <- ifelse(sif.df[, 4] < 0, "#56B4E9", "#CC79A7")
    nNodes <- length(igraph::V(netTmp))
   
    layoutTmp <- igraph::layout_(netTmp, with_fr()) 
    # plot(netTmp, edge.arrow.size=.1, edge.color=igraph::E(netTmp)$color,
    #      #vertex.color=MS2netColsSub, 
    #      vertex.label.font=1, vertex.label.color= "gray83", vertex.label='',
    #       vertex.size=4, #vertex.shape=vertexShapesSub,
    #      vertex.label.cex=0)#, xlim=xlimTmp, ylim=ylimTmp) #layout=layout

    message(nNodes,
            " nodes with ", length(igraph::E(netTmp)), " edges identified at a spectral similarity (dot product score) >= ", minDotProdThresh, '\nOf the edges:\n',
            sum(sif.df[, 4] > 0), ' were based on a shared ion fragment pattern\n',
            sum(sif.df[, 4] < 0), ' were based on a shared neutral loss pattern\n')
    flush.console()
    object@network$specSimGraph <- netTmp 
    object@network$specSimLayout <- layoutTmp
    return(object)
}) # end function