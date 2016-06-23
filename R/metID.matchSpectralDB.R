#' match spectra to spectral databases in the .msp file format (NIST)
#'  
#' @param object A "CompMS2" class object.  
#' @param mspFile character string to an online or local .msp spectrum database file.
#' (default = "http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/LipidBlast_Posi_Plasma_vs2.msp")
#' @param minDBDotProdThresh minimum dot product spectral similarity score (default = 0.8).
#' @param ppmMS1 numeric minimum mass accuracy (ppm) between database (.msp)
#'  precursor masses and the CompMS2 composite spectrum MS1 m/z. 
#'  Dot products will only be calculated for database entries within this mass 
#'  accuracy threshold (default=10ppm).
#' @param binSizeMS2 numeric bin size for between database (.msp) fragment masses and the CompMS2 composite spectrum fragment ions (default = 0.1). 
#' @param minPropEx numeric (0-1) minimum proportion of total signal explained by the reference spectrum. If the minimum dot product score is not satisfied then a spectra which is above the minimum proportion explained will be matched. This is incorporated as remaining background signals can affect the dot product score (default = 0.6).
#' @return "CompMS2" class object with any database matches above the minimum dot product score. 
#' 
#' @export
setGeneric("metID.matchSpectralDB", function(object, ...) standardGeneric("metID.matchSpectralDB"))

setMethod("metID.matchSpectralDB", signature = "CompMS2", function(object, mspFile='http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/LipidBlast_Posi_Plasma_vs2.msp', minDBDotProdThresh=0.8, ppmMS1=10, binSizeMS2=0.1, minPropEx=0.6){
    # error handling
    stopifnot(!is.null(object))
    if(class(object) != "CompMS2"){
      stop('argument object must be a "CompMS2" class object')
    }
  # add parameters to object
  object@Parameters$minDBDotProdThresh <- minDBDotProdThresh
  object@Parameters$ppmMS1 <- ppmMS1
  object@Parameters$binSizeMS2 <- binSizeMS2
  object@Parameters$minPropEx <- minPropEx
  # empty results list if necessary
  if(length(object@spectralDB) == 0){
  object@spectralDB <- vector('list', length(object@compSpectra))  
  names(object@spectralDB) <- names(object@compSpectra)
  } 
  message('Reading lines from .msp file...Please wait\n')
  flush.console()
  
  msp <- readLines(mspFile)
  # remove empty lines
  msp <- msp[msp != '']
  nEntries <- grep('^NAME:', msp, ignore.case = T)
  message(length(nEntries), ' entries in .msp file\n')
  flush.console()
  
  splitFactorTmp <- rep(1:length(nEntries), diff(c(nEntries, length(msp) + 1)))
  # matrix of masses and intensities
  massIntIndx <- grep('^[0-9]', msp)
  massesInts <- unlist(strsplit(msp[massIntIndx], '\t| '))
  massesInts <- as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
  massesTmp <-  massesInts[seq(1, length(massesInts), 2)]
  relIntTmp <-  massesInts[seq(2, length(massesInts), 2)]
  # label with entry numbers
  names(massesTmp) <- splitFactorTmp[massIntIndx]
  names(relIntTmp) <- splitFactorTmp[massIntIndx]
  # entry info
  entryInfoIndx <- setdiff(1:length(msp), massIntIndx)
  entryInfo <- msp[entryInfoIndx]
  names(entryInfo) <- paste0(splitFactorTmp[entryInfoIndx], '_', gsub(':.+', '', entryInfo))
  entryInfo <- gsub('^[A-Z]*: ', '', entryInfo, ignore.case = T)
  DBprecursorMasses <- as.numeric(entryInfo[grep('PRECURSORMZ', names(entryInfo), ignore.case = T)])
  # add database data to table
  DBdf <- cbind(mass=massesTmp, Rel_Intensity=relIntTmp)
  dbNamesTmp <- paste0('DB_', row.names(DBdf))
  message('\ncalculating spectral similarities (dot product >= ', round(minDBDotProdThresh, 1),
          ' or minimum proportion of spectrum explained >=', minPropEx, ') between database and composite spectra...\n')
  flush.console()
  # all constituent spectra
  pb <- txtProgressBar(max=length(object@compSpectra), style = 3)
  noMatchesV <- as.numeric()
  for(i in 1:length(object@compSpectra)){
    setTxtProgressBar(pb, i)
    massIntTmp <- object@compSpectra[[i]][, c('mass', 'Rel_Intensity')]
    ms1MzIndx <- grep('MS1_mz', names(object@metaData[[i]]), ignore.case=T)
    ms1MzTmp <- unlist(object@metaData[[i]][ms1MzIndx])[1]
    dbMatchMasses <- which(abs(((ms1MzTmp - DBprecursorMasses)/ ms1MzTmp) * 1E06) <= ppmMS1)
    if(length(dbMatchMasses) == 0){
    noMatchesV <- c(noMatchesV, i)  
    next
    } 
    indxTmp <- row.names(DBdf) %in% dbMatchMasses 
    DBdfTmp <- DBdf[indxTmp, , drop=F]
   
    spectrumDBTmp <- rbind(massIntTmp, DBdfTmp)
    specNamesVecTmp <- c(rep('DB_0', nrow(massIntTmp)), dbNamesTmp[indxTmp])
    # padded integer labels
    
    labelsTmp <- paste0('(', seq(binSizeMS2, (1000 - binSizeMS2), binSizeMS2), ',', seq((2 * binSizeMS2), 1000, binSizeMS2), ']')
    massBinsIndivTmp <- cut(spectrumDBTmp[, 1], breaks=seq(binSizeMS2, 1000, binSizeMS2), labels=labelsTmp)   
    # empty bins
    indivSpecVec <- tapply(spectrumDBTmp[, 2], paste0(specNamesVecTmp, massBinsIndivTmp), sum)
    # identify any absent bins
    allBinNames <- paste0(rep(unique(specNamesVecTmp), each=length(labelsTmp)), rep(labelsTmp, length(unique(specNamesVecTmp))))
    # add absent bins as zeros
    allBinsTmp <- rep(0, length(allBinNames))
    names(allBinsTmp) <- allBinNames
    # ensure indivSpecVec is in right order
    allBinsTmp[match(names(indivSpecVec), allBinNames)] <- indivSpecVec
    
    indivSpecMat <- matrix(allBinsTmp, byrow=F, nrow=length(labelsTmp))
    # proportion explained
    propTicEx <- apply(indivSpecMat, 2, function(x){ 
      tmpIndx <- x > 0 
      compSpecIndxTmp <- indivSpecMat[, 1] > 0
      return(sum(indivSpecMat[tmpIndx & compSpecIndxTmp, 1])/ sum(indivSpecMat[compSpecIndxTmp, 1]))})
    # mean all pairwise dotproducts
    # dotProdMat <- t(indivSpecMat) %*% indivSpecMat
    dotProdMat <- crossprod(indivSpecMat)
    sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), nrow=nrow(dotProdMat), 
                            ncol=ncol(dotProdMat), byrow = T) 
    
    dotProdsTmp <- dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
    # first column is the spectrum dot prods
    matchIndxTmp <- dotProdsTmp[, 1] >= minDBDotProdThresh | propTicEx >= minPropEx
    dotProdsTmp <- dotProdsTmp[matchIndxTmp, , drop=F]
    propTicEx <- propTicEx[matchIndxTmp]
    dotProdsTmp <- dotProdsTmp[-1, 1]
    propTicEx <- propTicEx[-1]
    
    if(length(dotProdsTmp) > 0){
      DBmatchesTmp <- dbMatchMasses[matchIndxTmp[-1]]
      dupIndx <- duplicated(DBmatchesTmp) == F
      DBmatchesTmp <- DBmatchesTmp[dupIndx]
      dotProdsTmp <- dotProdsTmp[dupIndx]
      indxTmp <- rownames(DBdf) %in% DBmatchesTmp
      dbMassesInt <- DBdf[indxTmp, , drop=F]
      indxTmp <- gsub('_.+', '', names(entryInfo)) %in% DBmatchesTmp
      entryInfoTmp <- split(entryInfo[indxTmp], gsub('_.+', '', names(entryInfo[indxTmp])))
      compound_msp <- entryInfo[indxTmp][grep('_NAME', names(entryInfo[indxTmp]), ignore.case = T)]
      compound_msp <- paste0(compound_msp, '__',
                             gsub('\\.msp$', '', basename(mspFile)))
      names(entryInfoTmp) <- compound_msp
      compound_msp <- unlist(mapply(rep, compound_msp, each=table(row.names(dbMassesInt)), SIMPLIFY = F))
      
      dotProdsTmp <- unlist(mapply(rep, dotProdsTmp, each=table(row.names(dbMassesInt)), SIMPLIFY = F))
      propTicEx <- unlist(mapply(rep, propTicEx, each=table(row.names(dbMassesInt)), SIMPLIFY = F))
      dbMassesInt <- cbind(dbMassesInt, compound_msp, dotProductScore=dotProdsTmp, propTIC_explained=propTicEx)
      if(length(object@spectralDB[[i]]) > 0){
      object@spectralDB[[i]]$dbSpectra <- rbind(object@spectralDB[[i]]$dbSpectra, 
                                                dbMassesInt)
      object@spectralDB[[i]]$entryInfo <- c(object@spectralDB[[i]]$entryInfo, 
                                            entryInfoTmp)
      } else {
      object@spectralDB[[i]] <- list(dbSpectra=dbMassesInt, entryInfo=entryInfoTmp)
      }
      # message(i, '\n', max(dotProdsTmp[2:nrow(dotProdsTmp), 1]), '\n')
    } else {
      noMatchesV <- c(noMatchesV, i)
    }
  }
  # set the difference between those with any previous matches and new no matches vector
  noMatchesV <- setdiff(noMatchesV, which(sapply(object@spectralDB, length) > 0))
  
  message('\n', length(object@compSpectra)-length(noMatchesV), ' composite spectra (', 
          round(((length(object@compSpectra)-length(noMatchesV))/length(object@compSpectra)) * 100, 1), 
          '%) currently matched to spectral databases\n')
  flush.console()
  return(object)
}) # end function