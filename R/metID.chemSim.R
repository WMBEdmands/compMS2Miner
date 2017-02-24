#' tanimoto chemical similarity first network neighbours
#' 
#' @param object a "compMS2" class object.
#' @param autoPossId logical if TRUE the function will automatically add the name
#'  of the top annotation based on mean maximum 1st neighbour chemical similarity
#'  above the minimum chemical similarity score (default = FALSE). Caution if TRUE
#'  this will overwrite any existing possible_identities in the "metID comments"
#'  table. This functionality is intended as an automatic annotation identification tool prior to thorough examination of the data in \code{\link{compMS2Explorer}}.
#'  The intention is that automatic annotations can be used in the metID.rtPred
#'  retention prediction function as part of a seamless first-pass workflow.
#'  
#' @param minSimScore numeric must be values between 0-1 minimum tanimoto chemical similarity score (default = 0.8). Any mean maximum 1st neighbour chemical 
#' similarity scores will be considered for automatic possible identity 
#' addition to the "metID comments"
#'  table, the annotation with the highest mean maximum 1st neighbour chemical similarity score will then be automatically added to the "metID comments" table
#'  in the \code{\link{compMS2Explorer}}.
#'  
#' @param possContam numeric how many times does a possible annotation have
#'  to appear in the automatically generated possible annotations for it to be
#'  considered a contaminant and therefore not added to the "metID comment" table (default = 3, i.e. if a database name appears more than 3 times in the 
#'  automatic annotation table it will be removed).
#'  
#' @param bitsChemFP numeric values between 1024-4096 number of most frequent
#' atom-pairs in the DrugBank database see the ChemmineR function \code{\link{desc2fp}} for more details (default = 1024).
#' @param minEdges numeric minimum number of edges (i.e. connected/adjacent nodes/spectra) 
#' to consider a node/spectrum for chemical similarity (default=2). 
#' This filtration is performed after removal of isobaric spectra.
#' For example nodes with only one edge/adjacent node are more likely to produce 
#' false-positive annotations. For more robust chemical similarity scoring 
#' consider increasing this number.
#' @param verbose logical if TRUE display progress bars.
#'  
#' @details this function can only be utilized after running \code{\link{metID.corrNetwork}} and/or \code{\link{metID.specSimNetwork}}. The purpose of this function is to provide first-pass automatic metabolite annotation. The tanimoto chemical similarity score is first calculated from a 1024-4096 bit chemical fingerprint for every best annotation SMILES code. For annotations of each composite spectrum the maximum chemical similarity score with any first neighbours (either by correlation and/or spectral similarity) are identified and the weighted arithmetic mean maximum chemical similarity score of 1st neighbours calculated. The mean is weighted based on the mean spectral similarity and/or
#' correlation coefficient value for an edge pair (i.e. if two spectra are connected by both spectral similarity and
#' correlation then a mean value of the two will be calculated and used as the weight).
#' This is to ensure that the more similar or highly correlated two composite spectra
#' are the higher the contribution to the maximum chemical similarity score. 
#' A new column "MMNNCSS" is added to the best annotation tables for any composite spectra with at least one composite spectrum network neighbour. Furthermore, the best Annotation table is
#' sorted according to this new likely annotation score, this give the user a rapid
#' means to establish a likely annotation based on chemical similarity with neighbouring node annotations.
#'  
#'  Optionally, the top annotations for a composite spectrum can be automatically 
#'  added to the "metID comments" table in the compMS2Explorer application. If the 
#'  argument \strong{autoPossId} is TRUE (default = FALSE) the function will automatically add the name
#'  of the top annotation based on mean maximum 1st neighbour chemical similarity
#'  above the minimum chemical similarity score (argument minSimScore, default = 0.8).
#' @examples 
#' library(compMS2Miner)
#' compMS2Example <- metID(compMS2Example, 'chemSim', minSimScore=0.8, 
#'                         autoPossId=TRUE)
#' @export 
setGeneric("metID.chemSim", function(object, ...) standardGeneric("metID.chemSim"))

setMethod("metID.chemSim", signature = "compMS2", function(object, autoPossId=FALSE,
                                                           minSimScore=0.8,
                                                           possContam=3, 
                                                           bitsChemFP=1024,
                                                           minEdges=2, 
                                                           verbose=TRUE, ...){
  # error handling
  if(!require(ChemmineR)){
    stop('package ChemmineR must be installed to use this function')
  }
  if(!require(ChemmineOB)){
    stop('package ChemmineOB must be installed to use this function')
  }
  if(!require(fingerprint)){
    stop('package fingerprint must be installed to use this function')
  }
  if(!require(igraph)){
    stop('package igraph must be installed to use this function')
  }
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  if(length(DBanno(object)) == 0){
    stop('The function metID.dbAnnotate and metID.dbProb must be run before chemical similarity calculation')
  } 
  if(length(BestAnno(object)) == 0){
    stop('The function metID.dbProb must be run before chemical similarity calculation')  
  }
  
  if(length(network(object)) == 0){
    stop('either the correlation network and/or spectral similarity network must be calculated before chemical similarity calculation')
  }
  # add parameters
  Parameters(object)$autoPossId_chemSim <- autoPossId
  Parameters(object)$minSimScore <- minSimScore
  Parameters(object)$possContam <- possContam
  Parameters(object)$bitsChemFP <- bitsChemFP
  
  # index which bestAnnos have entries
  emptyIndx <- sapply(BestAnno(object), nrow) > 0
  if(all(emptyIndx == FALSE)){
    stop('none of composite spectra have any dbAnnotate matches')
  }
  # match bestAnno names to vertices
  netGraphIndx <- names(network(object)) %in% c("corrNetworkGraph", "specSimGraph")
  vertNamesTmp <- unlist(lapply(network(object)[netGraphIndx], function(x) igraph::V(x)$name))
  vertNamesTmp <- unique(vertNamesTmp)
  netBestIndx <- match(names(BestAnno(object)), vertNamesTmp)
  # edge pairs from network
  edgePairs <- unlist(lapply(network(object)[netGraphIndx], function(x) attr(igraph::E(x), 'vnames')))
  # edge values i.e. corr coeff or spectral similarities for weighting
  edgeValues <- unlist(lapply(network(object)[netGraphIndx], function(x) abs(igraph::E(x)$value)))
  # duplicate edges 
  meanEdgeValues <- tapply(edgeValues, edgePairs, mean)
  edgePairs <- unique(edgePairs)
  # sort meanEdgeValues
  meanEdgeValues <- meanEdgeValues[match(edgePairs, names(meanEdgeValues))]
    
  vertIndxTmp <- netBestIndx[!is.na(netBestIndx) & emptyIndx]
  vertNamesSub <- vertNamesTmp[vertIndxTmp]
  # create vertex pairs list
  vertPairs <- lapply(vertNamesSub, function(x){
    # create regex to id correct edgepairs
    regExTmp <- paste0('^', x, '\\|', '|', '\\|', x, '$')
    edgePairsTmp <- edgePairs[grep(regExTmp, edgePairs)]
    edgePairValsTmp <- meanEdgeValues[grep(regExTmp, edgePairs)] 
    edgePairsTmp <- gsub(regExTmp, '', edgePairsTmp)
    dupIndx <- duplicated(edgePairsTmp) == FALSE
    edgePairValsTmp <- edgePairValsTmp[dupIndx]
    edgePairsTmp <- edgePairsTmp[dupIndx]
    indxTmp <- edgePairsTmp %in% vertNamesSub
    edgePairsTmp <- edgePairsTmp[indxTmp]
    edgePairValsTmp <- edgePairValsTmp[indxTmp]
    if(length(edgePairsTmp) > 0){
    edgePairsTmp <- c(x, edgePairsTmp)
    names(edgePairsTmp) <- c(0, edgePairValsTmp)
    return(edgePairsTmp)
    }})
  # null index i.e. composite spectrum connected to a composite spectrum with no
  # best annotations
  # vertNamesSub <- vertNamesSub[sapply(vertPairs, is.null) == FALSE]
  vertPairs <- vertPairs[sapply(vertPairs, is.null) == FALSE]
  # if too close in mass then remove
  mzPrec <- sapply(metaData(object)[emptyIndx], function(x) x[grep('_MS1_mz', names(x))][[1]][1])
  vertPairs <- lapply(vertPairs, function(x){
  remVertsTmp <- NULL 
  # greater than ppm than the 1st
  # tmpIndx <- abs({{mzPrec[x[1]] - mzPrec[x[-1]]}/mzPrec[x[1]]} * 1E06) > Parameters(object)$precursorPpm
  tmpIndx <- abs(mzPrec[x[1]] - mzPrec[x[-1]]) > 1.5
  if(sum(tmpIndx) > minEdges){
  remVertsTmp <- c(x[1], x[-1][tmpIndx])
  }
  return(remVertsTmp)
  })
  vertPairs <- vertPairs[sapply(vertPairs, is.null) == FALSE]
  # extract smiles codes and replace rownames to compSpecName
  smilesDf <- do.call(rbind, BestAnno(object)[emptyIndx])
  nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
  smilesDf$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[emptyIndx]),
                                    each=nrowBestAnno))
  
  # check if any commented then only retain these
  # alreadyAnno <- grepl('metID\\.matchSpectralDB', Comments(object)$user_comments)
  # keepIndx <- rep(TRUE, nrow(smilesDf))
  # if(any(alreadyAnno)){
  # indxTmp <- smilesDf$specNamesTmp %in%  Comments(object)$compSpectrum[alreadyAnno]
  # indxTmp <- indxTmp & {{smilesDf$DBname %in% Comments(object)$possible_identity[alreadyAnno]} == FALSE}
  # keepIndx[indxTmp] <- FALSE
  # smilesDf <- smilesDf[keepIndx, , drop=FALSE]
  # }
  # remove substr types
  smilesDf[is.na(smilesDf)] <- ''
  smilesDf <- smilesDf[smilesDf$SubStr_type == '', , drop=FALSE]
  # subset
  smilesDf <- smilesDf[, colnames(smilesDf) %in% c('DBid', 'SMILES', 'DBname',
                                                   'SubStr_type', 'ESI_type', 'WebAddress',
                                                   'chemFP', 'specNamesTmp')]
  # remove duplicates 
  smilesDfSub <- smilesDf[duplicated(smilesDf$DBid) == FALSE, , drop=FALSE]

if(('chemFP' %in% colnames(smilesDfSub)) == FALSE){
  # add new column
  smilesDfSub$chemFP <- ''
  Parameters(object)$bitsChemFP <- bitsChemFP
  # if possible and bitsChemFP equal to 1024 match as many ids as possible to internal
  # databases
  if(Parameters(object)$bitsChemFP == 1024){
    # add precalculated chemical fingerprints if available
    availDbs <- c('HMDB', 'LMSD', 'drugBank', 'T3DB', 'ReSpect')
    # attach 
    allDbIds <- do.call(c, lapply(availDbs, function(x) eval(parse(text=paste0('data(', x, ')')))))
    # extract dbids
    allDbIds <- do.call(c, lapply(availDbs, function(x) eval(parse(text=paste0(x, '$Unique_DB_ID')))))
    # extract chemical fingerprint
    names(allDbIds) <- do.call(c, lapply(availDbs, function(x){
      eval(parse(text=paste0(x, '$chemFP')))}))
    
    # add any matching database ids to smiles dataframe
    indxTmp <- match(smilesDfSub$DBid, allDbIds)
    smilesDfSub$chemFP[!is.na(indxTmp)] <- names(allDbIds)[indxTmp[!is.na(indxTmp)]]
    }
  # empty chemFP
  indxTmp <- smilesDfSub$chemFP == ''
  if(any(indxTmp)){  
  # convert first to SMILES set object and then SDF (ChemmineR/OB)
  message('Converting ', sum(indxTmp), ' SMILES codes -> SDF -> atom-pair descriptors -> ', bitsChemFP, 
          ' most common atom pairs (DrugBank.ca) chemical fingerprints this only needs to be performed once...please wait\n')
  flush.console()
  SDFtmp <- tryCatch({
    suppressWarnings(ChemmineR::smiles2sdf(gsub('\\*', '', smilesDfSub$SMILES[indxTmp])))
  }, error=function(cond) {
    return(0)
  }, warning=function(cond){
    return(0)})
  if(!is.numeric(SDFtmp)){
  SDFtmp@ID <- as.character(smilesDfSub$DBid[indxTmp])
  
  # if necessary alert user to failed apset converts
  errorIndx <- ChemmineR::validSDF(SDFtmp)
  
  if(any(errorIndx == FALSE)){
  message('The following database entries could not be converted to atom-pair ',
          'descriptors:\n', paste0(smilesDfSub$DBid[indxTmp][errorIndx == FALSE], ' ', paste0('(', smilesDfSub$DBname[indxTmp][errorIndx == FALSE], ')'), sep='\n'),
          'Chemical similarity scores cannot be calculated.')
    flush.console()
    # remove invalid sdf
    SDFtmp <- SDFtmp[errorIndx]
    # remove from smilesDfSub
    smilesDf <- smilesDf[-which(smilesDf$DBid %in% smilesDfSub$DBid[indxTmp][errorIndx == FALSE]), , drop=FALSE]
    smilesDfSub <- smilesDfSub[-which(indxTmp)[errorIndx == FALSE], , drop=FALSE]
    indxTmp <- smilesDfSub$chemFP == ''
  }
  # create APset atom pair descriptors files
  #   tryCatch({
  apsetTmp <- suppressWarnings(ChemmineR::sdf2ap(SDFtmp))
  
  # create n bit finger prints of apset using ChemmineR 
  if(length(apsetTmp) > 0){
  fpsetTmp  <-  ChemmineR::desc2fp(apsetTmp,  descnames=Parameters(object)$bitsChemFP,
                                   type="matrix")
  
  smilesDfSub$chemFP[indxTmp] <- apply(fpsetTmp, 1, function(x) 
                                       paste0(which(x == 1), collapse=';'))
  }
  # calculate tanimoto similarity
  # pmt <- proc.time()
  } else {
    message('None of the missing database entries could not be converted to atom-pair ',
            'descriptors:\n',
            'Chemical similarity scores cannot be calculated.')
    flush.console()
    errorIndx <- rep(FALSE, length(smilesDfSub$SMILES[indxTmp]))
    # remove from smilesDfSub
    smilesDf <- smilesDf[-which(smilesDf$DBid %in% smilesDfSub$DBid[indxTmp][errorIndx == FALSE]), , drop=FALSE]
    smilesDfSub <- smilesDfSub[-which(indxTmp)[errorIndx == FALSE], , drop=FALSE]
    indxTmp <- smilesDfSub$chemFP == '' 
  }
  } 
}
fpsetTmp <-  do.call(rbind, lapply(smilesDfSub$chemFP, function(x){
                    indxTmp <- as.numeric(strsplit(x, ';')[[1]])
                    return(ifelse(1:Parameters(object)$bitsChemFP %in% indxTmp, 1, 0))}))
row.names(fpsetTmp) <- smilesDfSub$DBid
message('calculating tanimoto similarity scores ', prettyNum(nrow(smilesDfSub), big.mark = ','), 
        ' smiles...please wait.\n')
flush.console()
# split matrix in two equal chunks of 2000



################################################################################
# # from fingerprint package function .tanimoto.sim.mat 
# mat <- fpsetTmp %*% t(fpsetTmp)
# len <- length(fpsetTmp[, 1])
# s <- mat.or.vec(len, len)
# smilesTanimoto <- .C("m_tanimoto", as.double(mat), as.integer(len), as.double(s),
#                      PACKAGE = "fingerprint")
# smilesTanimoto <- matrix(smilesTanimoto[[3]], nrow=len, ncol=len, byrow=TRUE
################################################################################

# expand smiles Tanimoto to size of smilesDf
dbIdIndxTmp <- match(smilesDf$DBid, smilesDfSub$DBid)
# add chemical fingerprints to smilesDf
smilesDf$chemFP <- smilesDfSub$chemFP[dbIdIndxTmp] 

# smilesTanimoto <- smilesTanimoto[indxTmp, indxTmp]
# fpsetTmp <- fpsetTmp[indxTmp, ]
# specNamesTmp <- gsub('^CC_|\\.[^\\.]*$', '', row.names(smilesDf))
# colnames(smilesTanimoto) <- gsub('^CC_|\\.[^\\.]*$', '', row.names(smilesDf))
# row.names(smilesTanimoto) <- colnames(smilesTanimoto)
# calculate the mean chemical similarity of first neighbours parallel if necessary
# if(Parameters(object)$nCores > 0){
#     if(!require(foreach)){
#       stop('the foreach package must be installed to run a parallel process.\n')
#     }
#       pmt <- proc.time()
#       nCores <- Parameters(object)$nCores
#       actualCores <- parallel::detectCores()
#       nCores <- ifelse(nCores >= actualCores, actualCores, nCores)
#       message(paste0("Starting SNOW cluster with ", nCores, " local sockets..."))
#       flush.console()
#       cl <- parallel::makeCluster(nCores)
#       doSNOW::registerDoSNOW(cl)
#       specNamesTmp <- smilesDf$specNamesTmp
#       # foreach and dopar from foreach package
#       chemSim1stNeigh <- foreach(x=vertPairs, .packages=c('fingerprint')) %dopar% {# from fingerprint package function .tanimoto.sim.mat 
#           indxTmp <- match(specNamesTmp, x)
#           dbIdIndxTmpSub <- dbIdIndxTmp[!is.na(indxTmp)]
#           fpsetTmpSub <- fpsetTmp[unique(dbIdIndxTmpSub), , drop=FALSE]
#           # mat <- fpsetTmpSub %*% t(fpsetTmpSub)
#           mat <- tcrossprod(fpsetTmpSub)
#           len <- length(fpsetTmpSub[, 1])
#           s <- mat.or.vec(len, len)
#           smilesTanimoto <- .C("m_tanimoto", as.double(mat), as.integer(len), 
#                                as.double(s), PACKAGE = "fingerprint")
#           smilesTanimoto <- matrix(smilesTanimoto[[3]], nrow=len, ncol=len, 
#                                    byrow=TRUE)
#           namesTmp <-  specNamesTmp[!is.na(indxTmp)]
#           expandIndxVert <- as.numeric(as.factor(dbIdIndxTmpSub))
#           expandIndxNeigh <- expandIndxVert[namesTmp %in% x[2:length(x)]]
#           expandIndxVert <- expandIndxVert[namesTmp %in% x[1]]
#           smilesTanimoto <- smilesTanimoto[expandIndxNeigh, expandIndxVert, drop=F]
#           row.names(smilesTanimoto) <- namesTmp[namesTmp %in% x[2:length(x)]]
#           colnames(smilesTanimoto) <- namesTmp[namesTmp %in% x[1]]
#           # vector weights
#           weightsV <- as.numeric(names(x))[unique(match(row.names(smilesTanimoto), x))] 
#           # smilesTanimoto <- smilesTanimoto[row.names(smilesTanimoto) %in% x[2:length(x)], colnames(smilesTanimoto) %in% x[1], drop=FALSE]
#           smilesTanimoto[is.na(smilesTanimoto)] <- 0
#           meanChemSimTmp <- apply(smilesTanimoto, 2, function(y) weighted.mean(tapply(y, row.names(smilesTanimoto), max), w = weightsV))
#       }
#       # stop SNOW cluster
#       parallel::stopCluster(cl)
#       proc.time() - pmt
#       chemSim1stNeigh <- do.call(c, chemSim1stNeigh)
#       } else {
      # for loop
      chemSim1stNeigh <- vector('list', length(vertPairs))
      if(verbose == TRUE){ pb <- txtProgressBar(max=length(vertPairs), style=3)}  
      pmt <- proc.time()
      for(x in 1:length(vertPairs)){
        if(verbose == TRUE){setTxtProgressBar(pb, x)}
        vpTmp <- vertPairs[[x]]
        indxTmp <- match(smilesDf$specNamesTmp, vpTmp)
        # match order of smiles specNames and vertix pairs
        sortIndxTmp <- which(!is.na(indxTmp))
        sortIndxTmp <- sortIndxTmp[order(indxTmp[!is.na(indxTmp)])]
        # subset vpTmp
        vpIndx <- sort(unique(indxTmp[!is.na(indxTmp)]))
        vpTmp <- vpTmp[vpIndx]
        if(length(vpTmp) > 1){
        firstV <- which(names(vpTmp) == '0')
        vpTmp <- c(vpTmp[firstV], vpTmp[-firstV])
        # dbIdIndxTmpSub <- dbIdIndxTmp[!is.na(indxTmp)]
        # fpsetTmpSub <- fpsetTmp[unique(dbIdIndxTmpSub), , drop=FALSE]
        dbIdsTmp <- smilesDf$DBid[sortIndxTmp]
        names(dbIdsTmp) <- smilesDf$specNamesTmp[sortIndxTmp]
        fpsetTmpSub <- fpsetTmp[as.character(unique(dbIdsTmp)), , drop=FALSE]
        mat <- tcrossprod(fpsetTmpSub)
        len <- length(fpsetTmpSub[, 1])
        s <- mat.or.vec(len, len)
        smilesTanimoto <- .C("m_tanimoto", as.double(mat), as.integer(len), 
                             as.double(s), PACKAGE = "fingerprint")
        smilesTanimoto <- matrix(smilesTanimoto[[3]], nrow=len, ncol=len, 
                                 byrow=TRUE)
        colnames(smilesTanimoto) <- unique(dbIdsTmp)
        row.names(smilesTanimoto) <- unique(dbIdsTmp)
        # DBNameTmp <- smilesDf$DBid[!is.na(indxTmp)]
        # expandIndxVert <- as.numeric(as.factor(dbIdIndxTmpSub))
        expandIndxNeigh <- dbIdsTmp[names(dbIdsTmp) %in% vpTmp[2:length(vpTmp)]]
        expandIndxVert <- dbIdsTmp[names(dbIdsTmp) %in% vpTmp[1]]
        smilesTanimoto <- smilesTanimoto[as.character(expandIndxNeigh), as.character(expandIndxVert), drop=FALSE]
        # vector weights
        weightsV <- as.numeric(names(vpTmp))[-1] 
        # smilesTanimoto <- smilesTanimoto[row.names(smilesTanimoto) %in% vpTmp[2:length(vpTmp)], colnames(smilesTanimoto) %in% vpTmp[1], drop=FALSE]
        smilesTanimoto[is.na(smilesTanimoto)] <- 0
        chemSimTmp <- apply(smilesTanimoto, 2, function(y) weighted.mean(tapply(y, names(expandIndxNeigh), max), w = weightsV))
        
        if(length(chemSimTmp) > 0){
        names(chemSimTmp) <- paste0(names(expandIndxVert), ';', names(chemSimTmp))
        }
        chemSim1stNeigh[[x]] <- chemSimTmp
        }
      }  
      proc.time() - pmt
      chemSim1stNeigh <- do.call(c, chemSim1stNeigh)  
      # chemSim1stNeigh <- do.call(c, lapply(vertPairs, function(x){
      #   # from fingerprint package function .tanimoto.sim.mat 
      #   indxTmp <- match( smilesDf$specNamesTmp), x)
      #   fpsetTmpSub <- fpsetTmp[!is.na(indxTmp), , drop=FALSE]
      #   mat <- fpsetTmpSub %*% t(fpsetTmpSub)
      #   len <- length(fpsetTmpSub[, 1])
      #   s <- mat.or.vec(len, len)
      #   smilesTanimoto <- .C("m_tanimoto", as.double(mat), as.integer(len), as.double(s),
      #                        PACKAGE = "fingerprint")
      #   smilesTanimoto <- matrix(smilesTanimoto[[3]], nrow=len, ncol=len, byrow=TRUE)
      #   row.names(smilesTanimoto) <- smilesDf$specNamesTmp[!is.na(indxTmp)])
      #   colnames(smilesTanimoto) <- row.names(smilesTanimoto)
      #   smilesTanimoto <- smilesTanimoto[row.names(smilesTanimoto) %in% x[2:length(x)], colnames(smilesTanimoto) %in% x[1], drop=FALSE]
      #   meanChemSimTmp <- apply(smilesTanimoto, 2, function(y) mean(tapply(y, row.names(smilesTanimoto), max)))
      # }))
# }
# add chem sim score back to best annotations
indxTmp <- match(paste0(smilesDf$specNamesTmp, ';', smilesDf$DBid), names(chemSim1stNeigh))
smilesDf$MMNNCSS <- 0
# smilesDf <- smilesDf[indxTmp, , drop=FALSE]
smilesDf$MMNNCSS[!is.na(indxTmp)] <- chemSim1stNeigh[indxTmp[!is.na(indxTmp)]]
# split to list
bestAnnoTmp <- split(smilesDf, f=smilesDf$specNamesTmp)

DBanno(object)[names(bestAnnoTmp)] <- lapply(names(bestAnnoTmp), function(x){
  bestAnnoDfTmp <- DBanno(object)[[x]]
  bestAnnoDfTmp$MMNNCSS <- 0
  bestAnnoDfTmp$chemFP <- ''
  indxTmp <- match(bestAnnoDfTmp$DBid, bestAnnoTmp[[x]]$DBid)
  # remove spec names
  bestAnnoTmp[[x]]$specNamesTmp <- NULL
  bestAnnoDfTmp$MMNNCSS[!is.na(indxTmp)] <- bestAnnoTmp[[x]]$MMNNCSS[indxTmp[!is.na(indxTmp)]]
  bestAnnoDfTmp$chemFP[!is.na(indxTmp)] <- bestAnnoTmp[[x]]$chemFP[indxTmp[!is.na(indxTmp)]]
  # sort by chem sim
  bestAnnoDfTmp <- bestAnnoDfTmp[order(bestAnnoDfTmp$MMNNCSS, decreasing = TRUE), , drop=FALSE]
})

# add empty columns to bestAnno where necessary
noChemFPindx <- sapply(DBanno(object)[emptyIndx], function(x) 'chemFP' %in% colnames(x)) == FALSE

if(any(noChemFPindx)){
  tmpSpecNames <- names(DBanno(object)[emptyIndx][noChemFPindx])
  DBanno(object)[tmpSpecNames] <- lapply(tmpSpecNames, function(x){
    bestAnnoDfTmp <- DBanno(object)[[x]]
    bestAnnoDfTmp$MMNNCSS <- 0
    bestAnnoDfTmp$chemFP <- ''
    return(bestAnnoDfTmp)
  }) 
}

if(autoPossId ==  TRUE){
  # message('CAUTION: do you wish to automatically add the most probable database annotations to the comments table (column possible_identity)?\n type (Y/N) and press [enter] to continue (N.B. matchSpectralDB matches will not be overwritten):')
  # flush.console()
  # overWritePossId <- readline()
  # if(overWritePossId == 'Y'){
    smilesDf <- smilesDf[smilesDf$MMNNCSS >= minSimScore, , drop=FALSE]
    if(nrow(smilesDf) > 0){
    # met Id comments table
    metIDcomments <- Comments(object)
    
    smilesDf <- smilesDf[order(smilesDf$MMNNCSS, decreasing = TRUE), , drop=FALSE]
    # remove duplicate composite spectra
    smilesDf <- smilesDf[duplicated(smilesDf$specNamesTmp) == FALSE, , drop=FALSE]
    # id possible contaminants
    possContaminants <- table(smilesDf$DBname)
    if(any(possContaminants > possContam)){
    indxTmp <- possContaminants > possContam
    message('The following automatic possible annotations have been identified as possible contaminants (i.e. occuring more than ', possContam, ' times see possContam argument) and will be named "possible_contaminant" in the Comments table and flagged "metID.chemSim":\n', paste0(1:sum(indxTmp), '. ', names(possContaminants)[indxTmp], ' (occurs n=', possContaminants[indxTmp], ' times)\n'))
    flush.console()
    contamIndx <- (smilesDf$DBname %in% names(possContaminants[indxTmp]))
    contamSpec <- smilesDf$specNamesTmp[contamIndx]
    # add to comments
    indxTmp <- metIDcomments$compSpectrum %in% contamSpec
    metIDcomments$user_comments[indxTmp] <- 'metID.chemSim'
    metIDcomments$possible_identity[indxTmp] <- 'possible_contaminant'
    smilesDf <- smilesDf[contamIndx == FALSE, , drop=FALSE]
    }
    
    message(nrow(smilesDf), ' composite spectra out of ', nrow(metIDcomments), ' total (', round((nrow(smilesDf)/nrow(metIDcomments)) * 100, 2), '%) have been annotated based on a mean maximum 1st neighbour chemical similarity score (>= ', round(minSimScore, 2), '). First neighbours of the ',  
            paste0(ifelse(names(network(object))[netGraphIndx] == 'specSimGraph', 'spectral similarity', 'correlation'), collapse = ' and '), ' network(s) were utilized.\n\n')
    flush.console()
    message('These annotations will now be added to the "metID comments" table in compMS2Explorer.\n')
    flush.console()
    # add possible annotations to metIDcomments
    alreadyAnno <- Comments(object)$compSpectrum[grepl('metID\\.matchSpectralDB', Comments(object)$user_comments)]
    if(length(alreadyAnno) > 0){
    smilesDf <- smilesDf[{smilesDf$specNamesTmp %in% alreadyAnno} == FALSE, , drop=FALSE]
    }
    indxTmp <- match(metIDcomments$compSpectrum, smilesDf$specNamesTmp)
    metIDcomments$possible_identity[!is.na(indxTmp)] <- smilesDf$DBname[indxTmp[!is.na(indxTmp)]]
    metIDcomments$ESI_type[!is.na(indxTmp)] <- smilesDf$ESI_type[indxTmp[!is.na(indxTmp)]]
    metIDcomments$user_comments[!is.na(indxTmp)] <- 'metID.chemSim'
    Comments(object) <- metIDcomments
    } else {
    warning('no mean maximum 1st neighbour chemical similarity scores were above the minimum similarity score of ', round(minSimScore, 2), ' consider reducing this parameter.\n')  
    }
  } # end autoPossId
# proc.time() - pmt
 # message('calculating pca model...please wait.\n')
 # flush.console()
 # calculating PCA model
 # pcaModel <- pcaMethods::pca(smilesTanimoto, ...)
 # # return scores
 # pcaScoresTmp <- pcaModel@scores
 # smilesDfSub$SubStr_type <- NULL
 # pcaScoresTmp <- cbind(smilesDfSub, pcaScoresTmp) 
 # 
 # colnames(smilesTanimoto) <- smilesDfSub$DBid
 # row.names(smilesTanimoto) <- smilesDfSub$DBid
 # 
 # # replace upper tri with zero
 # smilesTanimoto[upper.tri(smilesTanimoto, diag=TRUE] <- 0
 # 
 # aboveSimScore <- which(smilesTanimoto >= minSimScore, arr.ind = TRUE)
 # 
 # # unique db ids above sim score
 # aboveSimDbIds <- cbind(smilesDfSub$DBid[aboveSimScore[, 1]], 
 #                        smilesDfSub$DBid[aboveSimScore[, 2]])
 # specIdsDbIds <- tapply(gsub('\\.[^\\.]*$', '', row.names(smilesDf)), smilesDf$DBid,
 #                        function(x) paste0(x, collapse=';'))
 # 
 # message(length(unique(c(smilesDfSub$DBid[aboveSimScore[, 1]], 
 #                             smilesDfSub$DBid[aboveSimScore[, 2]]))), ' database entries out of ', nrow(smilesDfSub), ' total have a chemical similarity score >= ', minSimScore, ' (Tanimoto)\n')
 # flush.console()
 # 
 # netTmp <- igraph::graph(as.vector(t(aboveSimDbIds[, c(1, 2)])))
 # igraph::E(netTmp)$tanimotoScore <- smilesTanimoto[aboveSimScore]
 #    
 # # add edge colours to igraph object based on fragment or neutral loss
 # nNodes <- length(igraph::V(netTmp))
 # 
 # layoutTmp <- igraph::layout_(netTmp, with_fr()) 
 # indxTmp <- match(igraph::V(netTmp)$name, smilesDfSub$DBid)
 # layoutTmp <- cbind(layoutTmp, smilesDfSub[indxTmp, ])
 # indxTmp <- match(layoutTmp$DBid, names(specIdsDbIds))
 # layoutTmp$compSpec <- specIdsDbIds[indxTmp]
 # message(nNodes,
 #         " nodes with ", length(igraph::E(netTmp)),
 #         " edges identified at a chemical similarity score (tanimoto) >= ", 
 #         minSimScore, '\n')
 # flush.console()
 # network(object)$chemSimGraph <- netTmp 
 # network(object)$chemSimLayout <- layoutTmp
 return(object)
}) # end function
