#' msp file from compMS2 class object
#' 
#' @details this function converts a compMS2 class object to an msp database
#' file inserting experimental details and taking annotation data from the 
#' comments table and best annotations. This creates a means to database successive
#' experiments and potentially concatenate to a laboratories existing spectral 
#' database. In this way previously annotated spectra from other experiments 
#' can be matched to new datasets. A unique splash code is also generated for
#' each composite spectrum and added to each database entry in the msp file.
#' install splash R from GitHub using devtools:
#'  devtools::install_github("berlinguyinca/spectra-hash", subdir="splashR")
#' @param object a "compMS2" class object
#' @param studyName character a study name tag for your dataset.
#' @param existingMsp character a full path or web address of an existing msp database
#' file to concatenate to. The new entries will appear underneath the last entry 
#' of the existing msp datafile file. If this argument is supplied then this file will be
#' overwritten.
#' @param outputDir character full path to a directory to write the msp file to.
#' Default is to take the current working directory obtained from \code{\link{getwd}}.
#' @param onlyCommented logical (default = FALSE) should only metabolites recorded
#' in the comments table be included or if true all of the composite spectra.
#' 
#' @return the msp file will be written to a file tagged with the study name and date in the output directory.
#' @export
setGeneric("metID.compMS2toMsp", function(object, studyName=NULL, existingMsp=NULL,
                                          outputDir=getwd(), onlyCommented=FALSE, 
                                          ...) standardGeneric("metID.compMS2toMsp"))

setMethod("metID.compMS2toMsp", signature = "compMS2", function(object, studyName,
                                                                existingMsp,
                                                                outputDir, onlyCommented,
                                                                ...){ 
  # error handling
  # include instrumental parameters from datafile
  if(!require(splashR)){
    stop('Please install splashR from GitHub using devtools: devtools::install_github("berlinguyinca/spectra-hash", subdir="splashR") this will be used to generated a unique splash code (spectral hash codes) for each spectrum.')
  }
  if(is.null(studyName)){
    stop('The studyName argument must be supplied to tag the entries in the msp file.')
  }
  studyNameEntry <- paste0('STUDYNAME: ', studyName)
  # index which bestAnnos have entries
  emptyIndx <- sapply(BestAnno(object), nrow) > 0
  if(all(emptyIndx == FALSE)){
    stop('none of composite spectra have any dbAnnotate matches')
  }
  
  # extract smiles codes 
  bestAnnoDf <- do.call(rbind, BestAnno(object)[emptyIndx])
  nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
  bestAnnoDf$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[emptyIndx]),
                                               each=nrowBestAnno))
  # match best anno to comment
  if(nrow(Comments(object)) == 0){
    stop('No comments have yet been made and therefore no possible_identity have been entered in the column.\n\nEither run the automatic chemical similarity metabolite identification (see ?metID.chemSim or metID.buildConsensus) or visualize with ?compMS2Explorer and add putative annotation names to the possible_identity column in the "met ID comments" tabbed panel') 
  }
  message('Extracting identified compounds from comments table and generating splash codes.\n')
  flush.console()
  # extract comments table
  metIDcomments <- Comments(object)
  # match possible_identity to best anno table
  indxTmp <- match(bestAnnoDf$specNamesTmp, metIDcomments$compSpectrum)
  bestAnnoDf$possible_identity <- metIDcomments$possible_identity[indxTmp]
  bestAnnoDf$compound_class <- metIDcomments$compound_class[indxTmp]
  bestAnnoDf$user_comments <- metIDcomments$user_comments[indxTmp]
  # which dbNames and possible_identity names match
  bestAnnoDf$nameMatch <- bestAnnoDf$DBname == bestAnnoDf$possible_identity
  # sort by bc score if present
  if(!is.null(bestAnnoDf$BC_meanScore)){
    bestAnnoDf <- bestAnnoDf[order(bestAnnoDf$BC_meanScore, decreasing = TRUE), , drop=FALSE]
  }
  
  if(onlyCommented==TRUE){
  namesTmp <- unique(bestAnnoDf$specNamesTmp[bestAnnoDf$nameMatch])  
  } else {
  namesTmp <- names(compSpectra(object))
  }
  # create entries
  entriesList <- do.call(c, lapply(namesTmp, function(x){
    # name 
    indxTmp <- which(bestAnnoDf$nameMatch & {bestAnnoDf$specNamesTmp %in% x})
    nameTmp <- paste0('NAME: ', ifelse(length(indxTmp) > 0, bestAnnoDf$possible_identity[indxTmp[1]], 'not_identified'), 
                      ' ', ifelse(length(indxTmp) > 0, bestAnnoDf$ESI_type[indxTmp[1]], ''))
    
    # top 5 annotations
    topAnnos <- 'TOPANNO: no_annotations'
    topAnnoSmi <- 'no_annotations'
    topAnnoEsi <- ''
    topIndxTmp <- which(bestAnnoDf$specNamesTmp %in% x)
    if(length(topIndxTmp) > 0){
      nTopAnno <- ifelse(length(topIndxTmp) > 5, 5, length(topIndxTmp))
      topAnnos <- paste0('TOPANNO', 1:nTopAnno, ': ',
                         paste(bestAnnoDf$DBname[topIndxTmp[1:nTopAnno]],
                               bestAnnoDf$ESI_type[topIndxTmp[1:nTopAnno]]))
      topAnnoSmi <- bestAnnoDf$SMILES[topIndxTmp[1]]
      topAnnoEsi <- bestAnnoDf$ESI_type[topIndxTmp[1]]
    }
    modeTmp <- paste0('MODE: ', ifelse(Parameters(object)$mode == 'neg', 'negative', 'positive'))
    smilesTmp <- paste0('SMILES: ', 
                        ifelse(length(indxTmp) > 0, bestAnnoDf$SMILES[indxTmp[1]],
                               topAnnoSmi))
    
    atomF <- ''
    # generate atomic formula from smiles
    if(smilesTmp != ''){
      # remove all numbers and punctuation from smiles
      eleComp <- gsub('SMILES: |\n|[[:punct:]]|[[:digit:]]', '', smilesTmp)
      # capitalize aromatic lower case letters
      for(arom in c('c', 'o', 's', 'n')){
        eleComp <- gsub(arom, toupper(arom), eleComp) 
      }
      # split into elements
      eleComp <- gsub('([[:upper:]])', ' \\1', eleComp)
      eleComp <- strsplit(eleComp, ' ')[[1]]
      eleComp <- eleComp[eleComp != '']
      atomF <- table(eleComp)
      atomF <- paste0(names(atomF), ifelse(atomF == 1, '', atomF), collapse = '')
    }
    atomF <- paste0('FORMULA: ', atomF)
    webAddTmp <- paste0('DBWEBLINK: ', ifelse(length(indxTmp) > 0, 
                                              paste0(bestAnnoDf$WebAddress[indxTmp[1]],
                                                     gsub('_.+', '', bestAnnoDf$DBid[indxTmp[1]]),
                                                     collapse=''), ''))
    # metaData
    metaDataTmp <- metaData(object)[[x]]
    precMz <- paste0('PRECURSORMZ: ', 
                    ifelse(length(indxTmp) > 0, bestAnnoDf$candidateMass[indxTmp[1]],
                           metaDataTmp[grep('_MS1_mz', names(metaDataTmp))][[1]][1]))
    ms1Rt <- paste0('RETENTIONTIME: ', metaDataTmp[grep('_MS1_RT', names(metaDataTmp))][[1]][1])
    ms1Eic <- paste0('EICNO: ', metaDataTmp[grep('_MS1_EICno', names(metaDataTmp))][[1]][1])
    ms1Adduct <- gsub(' .+', '', metaDataTmp[grep('_MS1_adduct', names(metaDataTmp))][[1]][1])
    adductTmp <- ifelse(length(indxTmp) > 0, bestAnnoDf$ESI_type[indxTmp[1]], topAnnoEsi)
    adductTmp <- ifelse(adductTmp == '', ifelse(Parameters(object)$mode == 'neg', '[M-H]-', '[M+H]+'), 
                        adductTmp)
    adductTmp <- paste0('ADDUCT: ', adductTmp)
    userComm <- paste0('USERCOMMENTS: ', ifelse(length(indxTmp) > 0, bestAnnoDf$user_comments[indxTmp[1]], ''))
    ms1Adduct <- paste0('MS1ADDUCT: ', ms1Adduct)
    compClass <- paste0('COMPOUND_CLASS: ', ifelse(length(indxTmp) > 0, bestAnnoDf$compound_class[indxTmp[1]], '') )
      
    # mass spectrum
    specTmp <- compSpectra(object)[[x]]
    numPeaks <- paste0('NUMPEAKS: ', nrow(specTmp))
    peaksTmp <- paste0(round(specTmp$mass, 5), ' ', round(specTmp$intensity, 2))
    # splash code
    splashCode <- paste0('SPLASHCODE: ', getSplash(specTmp[, c('mass', 'intensity')]))
  
    subEntry <- c(nameTmp, studyNameEntry, precMz, adductTmp, modeTmp, smilesTmp, atomF,
                  splashCode, ms1Rt, ms1Eic, ms1Adduct, compClass, webAddTmp, userComm, topAnnos,
                  numPeaks, peaksTmp, '') # add empty lines
  }))
  
  if(!is.null(existingMsp)){
    message('Reading existing msp file and concatenating new entries.\n')
    flush.console()
    existMsp <- readLines(existingMsp)
    entriesList <- c(existMsp, entriesList)
    mspFileName <- existingMsp
  } else {
    mspFileName <- paste0(outputDir, '/', studyName, '_', gsub('-', '_', Sys.Date()), '.msp')
  }
  message('Saving .msp file containing (', length(namesTmp), ' new entries) in the following output directory:\n', mspFileName, '\n')
  flush.console()
  # write to output directory
  writeLines(paste0(entriesList, collapse='\n'), con=mspFileName)
}) # end function
