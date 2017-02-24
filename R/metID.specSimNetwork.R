#' CompMS2 spectral similarity network generation
#' 
#' @description generates a dot product spectral similarity network from the MS/MS fragmentation data. The resulting spectral similarity network can then be viewed in compMS2Explorer. Utilizes the \code{\link{graph}} function of the \code{\link{igraph}} package.
#' 
#' @param object A "compMS2" class object.  
#' @param minDotProdThresh minimum dot product spectral similarity score. 
#' If no value is supplied the default is to estimate
#' an optimal cut-off value based on network attributes from a series of spectral similarity cut-off values. This is based upon the relationship between the 
#' number of nodes, edges, number of clusters and the network density (i.e. n actual edges/n potential edges) at each spectral similarity cut-off value. A plot showing
#' the result of this estimation will be generated.
#' @param binSizeMS2 numeric MS2 bin size for spectral similarity matching (default = 0.1)
#' @param minClustSize numeric minimum number of connected nodes that is cluster size (default =3). If a cluster
#' of nodes is less than this number then it will be removed.
#' @return "compMS2" class object with an additional network graph of any peakTable features above the correlation threshold.
#' @references Koh Aoki, Yoshiyuki Ogata, and Daisuke Shibata
#' Approaches for Extracting Practical Information from Gene Co-expression Networks in Plant Biology
#' Plant Cell Physiol (2007) 48 (3): 381-390 first published online January 23, 2007 doi:10.1093/pcp/pcm013
#' 
#' @export
setGeneric("metID.specSimNetwork", function(object, ...) standardGeneric("metID.specSimNetwork"))

setMethod("metID.specSimNetwork", signature = "compMS2", function(object, 
                                                                  minDotProdThresh=NULL,  
                                                                  binSizeMS2=0.1,
                                                                  minClustSize=3){
    # error handling
    stopifnot(!is.null(object))
    if(class(object) != "compMS2"){
      stop('argument object must be a "compMS2" class object')
    }
  if(!require(igraph)){
    stop('the package igraph must be installed to utilize this function.')
  }
    
  # all constituent spectra
  allSpecTmp <- do.call(rbind, lapply(compSpectra(object), function(x){
   massIntTmp <- x[, c('mass', 'intensity')]
   #  preFragIntTmp <- x[, c('Precursorfrag.diff', 'intensity')]
   #  preFragIntTmp[, 1] <- as.numeric(preFragIntTmp[, 1])
   # preFragIntTmp <- preFragIntTmp[preFragIntTmp[, 1] > 2, , drop=FALSE]
   # # colnames(preFragIntTmp)[1] <- 'mass'
   #  return(preFragIntTmp)
   }))
  nrowBestAnno <- sapply(compSpectra(object), nrow)
  specNamesVecTmp <- do.call(c, mapply(rep, names(compSpectra(object)),
                                       each=nrowBestAnno))
  message('Fragment ions:\n')
  flush.console()
   fragsDotProds <- dotProdMatrix(as.matrix(allSpecTmp), specNamesVecTmp, binSizeMS2=binSizeMS2)

    # 2. precursor - neutral losses
    
    # all constituent spectra
    allSpecTmp <- lapply(compSpectra(object), function(x){
      #massIntTmp <- x[, c('mass', 'intensity')]
      preFragIntTmp <- x[, c('Precursorfrag.diff', 'intensity'), drop=FALSE]
      preFragIntTmp[, 1] <- as.numeric(preFragIntTmp[, 1])
      preFragIntTmp <- preFragIntTmp[preFragIntTmp[, 1] > 2, , drop=FALSE]
      if(nrow(preFragIntTmp) == 0){
      preFragIntTmp <- data.frame(Precursorfrag.diff=2, intensity=1, stringsAsFactors = FALSE)  
      }
      # preFragIntTmp <- cbind(preFragIntTmp, name=names(compSpectra(object))[x])
      return(preFragIntTmp)
      # }
      # colnames(preFragIntTmp)[1] <- 'mass'
    })
    
    nrowBestAnno <- sapply(allSpecTmp, nrow)
    specNamesVecTmp <- do.call(c, mapply(rep, names(compSpectra(object)),
                                         each=nrowBestAnno))
    allSpecTmp <- do.call(rbind, allSpecTmp)
    
    message('Neutral losses:\n')
    flush.console()
    preFragDotProds <- dotProdMatrix(as.matrix(allSpecTmp), specNamesVecTmp, binSizeMS2=binSizeMS2)
    
    # replace upper tri with zero
    fragsDotProds[upper.tri(fragsDotProds, diag=TRUE)] <- 0
    # replace upper tri with zero
    preFragDotProds[upper.tri(preFragDotProds, diag=TRUE)] <- 0
    
    # which dot product is larger fragment similarity or neutral loss?
    indxTmp <- preFragDotProds > fragsDotProds 
    # replace values, remove negatives and add prefragdotprods necessary and make 
    # prefrag dot products negative to distinguish them.
    fragsDotProds[fragsDotProds < 0] <- 0
    fragsDotProds[indxTmp] <- -preFragDotProds[indxTmp]
    # identify optimum correlation cutoff
    if(is.null(minDotProdThresh)){
      cat('Estimating optimal spectral similarity cut-off based on network attributes from a sequence of cut-off values.\n')
      minDotProdThresh <- optimCutOff(fragsDotProds)$estCutOff
      cat('estimate of optimal minDotProdThresh value:', minDotProdThresh, '.\n')
    }
    # add parameters to object
    Parameters(object)$minDotProdThresh <- minDotProdThresh
    
    # make preFragDot
    # calculate pca model
    # pcDp <- pcaMethods::pca(preFragDotProds, method="svd", nPcs=2, cv='q2', center = FALSE)
    
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
    igraph::E(netTmp)$value <- sif.df[, 4]
    
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
