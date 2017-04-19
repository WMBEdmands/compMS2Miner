#' match spectra to spectral databases in the .msp file format (NIST)
#' 
#' @details If the msp file contains smiles codes then any matches above the min
#' dot prod or minimum proportion of the spectrum explained will be added to the 
#' best annotations table. 
#' @param object A "compMS2" class object.  
#' @param mspFile character string to an online or local .msp spectrum database file.
#' (default = NULL, dependendent on polarity either a positive or negative mode
#' version of massbank hosted on github will be downloaded) See \url{https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MassBank_MSMS_Pos_Rev173_vs1.msp} and \url{https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MassBank_MSMS_Neg_Rev173_vs1.msp}. Additional .msp files can be downloaded
#' from the massbank of North America website \url{http://mona.fiehnlab.ucdavis.edu/downloads}.
#' @param minDBDotProdThresh minimum dot product spectral similarity score (default = 0.65).
#' @param ppmMS1 numeric minimum mass accuracy (ppm) between database (.msp)
#'  precursor masses and the CompMS2 composite spectrum MS1 m/z. 
#'  Dot products will only be calculated for database entries within this mass 
#'  accuracy threshold (default=10ppm).
#' @param binSizeMS2 numeric bin size for between database (.msp) fragment masses and the CompMS2 composite spectrum fragment ions (default = 0.1). 
#' @param autoPossId logical if TRUE and if the .msp file database entries contain
#' at the very least a SMILES code (for downstream metabolite identification purposes) 
#' all msp file database entries will be added to the best annotations table.
#' and the function will also automatically add the name of the best annotation 
#' (highest dot product score above the user defined threshold) 
#' to the comments table and a note ("metID.matchSpectralDB") added to the 
#' comments column. This represents a form of automated metabolite identification
#' which can still be modified in the comments table upon visualization with
#' compMS2Explorer.
#' @param verbose logical if TRUE display progress bars.
#' @return "compMS2" class object with any database matches above the minimum dot product score. 
#' 
#' @export
setGeneric("metID.matchSpectralDB", function(object, ...) standardGeneric("metID.matchSpectralDB"))

setMethod("metID.matchSpectralDB", signature = "compMS2", function(object, mspFile=NULL, 
                                                                   minDBDotProdThresh=0.65, 
                                                                   ppmMS1=10, 
                                                                   binSizeMS2=0.1, 
                                                                   # minPropEx=0.8, 
                                                                   autoPossId = TRUE,
                                                                   verbose=TRUE){
  # error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  if(is.null(mspFile)){
  mspFile <- ifelse(Parameters(object)$mode == 'pos', 
                   'https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MassBank_MSMS_Pos_Rev173_vs1.msp',
                   'https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MassBank_MSMS_Neg_Rev173_vs1.msp')
  }
  # add parameters to object
  Parameters(object)$minDBDotProdThresh <- minDBDotProdThresh
  Parameters(object)$ppmMS1 <- ppmMS1
  Parameters(object)$binSizeMS2 <- binSizeMS2
  # Parameters(object)$minPropEx <- minPropEx
  # empty results list if necessary
  if(length(spectralDB(object)) == 0){
  spectralDB(object) <- vector('list', length(compSpectra(object)))  
  names(spectralDB(object)) <- names(compSpectra(object))
  } 
  message('Reading lines from .msp file...Please wait\n')
  flush.console()
  
  msp <- readLines(mspFile)
  # remove empty lines
  msp <- msp[msp != '']
  nEntries <- grep('^NAME:', msp, ignore.case = TRUE)
  message(prettyNum(length(nEntries), big.mark = ','), ' entries in .msp file\n')
  flush.console()
  
  splitFactorTmp <- rep(1:length(nEntries), diff(c(nEntries, length(msp) + 1)))
  # matrix of masses and intensities
  massIntIndx <- which(grepl('^[0-9]', msp) & !grepl(': ', msp))
  massesInts <- unlist(strsplit(msp[massIntIndx], '\t| '))
  # remove any containing 
  massesInts <- as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
  # if any NAs remove from indx
  massesTmp <-  massesInts[seq(1, length(massesInts), 2)]
  relIntTmp <-  massesInts[seq(2, length(massesInts), 2)]
  # label with entry numbers
  names(massesTmp) <- splitFactorTmp[massIntIndx]
  names(relIntTmp) <- splitFactorTmp[massIntIndx]
  # entry info
  entryInfoIndx <- setdiff(1:length(msp), massIntIndx)
  entryInfo <- msp[entryInfoIndx]
  names(entryInfo) <- paste0(splitFactorTmp[entryInfoIndx], '_', gsub(':.+', '', entryInfo))
  entryInfo <- gsub('.+: ', '', entryInfo, ignore.case = TRUE)
  # remove duplicate entries
  entryInfo <- entryInfo[duplicated(names(entryInfo)) == FALSE]
  
  # identify pertinent msp file entries
  mandEntryData <- data.frame(matrix('', nrow=length(nEntries), ncol=7), stringsAsFactors = FALSE)
  colnames(mandEntryData) <- c('DBname', 'DBIonMode', 'DBSMILES', 'DBINCHI', 
                               'DBprecursorMasses', 'DBesiType', 'chState')
  # name of entry
  mandFieldsRegEx <- c('_NAME$', '_ION MODE$|_MODE$', '_SMILES$', '_INCHI$', 
                       '_PRECURSORMZ$|_PRECURSOR M/Z$|_PRECURSOR MZ$|_PEPMASS$',
                       '_PRECURSORTYPE$|_PRECURSOR TYPE$|_ADDUCT$|_ION TYPE$|_IONTYPE$',
                       '_CHARGE$')
  names(mandFieldsRegEx) <- colnames(mandEntryData)
  for(mF in 1:length(mandFieldsRegEx)){
    entInfoTmp <- entryInfo[grep(mandFieldsRegEx[mF], names(entryInfo), ignore.case = TRUE)]
    entryNosTmp <- gsub(mandFieldsRegEx[mF], '', names(entInfoTmp), ignore.case=TRUE)
    # remove any double entries
    dupIdx <- duplicated(entryNosTmp) == FALSE
    entInfoTmp <- entInfoTmp[dupIdx]
    entryNosTmp <- entryNosTmp[dupIdx]
    indxTmp <- row.names(mandEntryData) %in% entryNosTmp
    mandEntryData[indxTmp, names(mandFieldsRegEx[mF])] <- entInfoTmp
  }
  
  # slashes in adduct esi type names
  mandEntryData$DBesiType <- gsub('/.+', '', mandEntryData$DBesiType)
  # 1. ionization mode
  
  mfIdxTmp <- mandEntryData$DBIonMode == ''
  if(any(mfIdxTmp)){
    # use adduct representation
    mandEntryData$DBIonMode[mfIdxTmp] <- ifelse(grepl('-$', mandEntryData$DBesiType[mfIdxTmp]), 'neg', 'pos')
  }
  # convert to logical
  # DBIonMode <- grepl(Parameters(object)$mode, DBIonMode, ignore.case = TRUE) 
  
 if(all(mandEntryData$DBIonMode ==  '')){
    stop('.msp file contains no ', ifelse(Parameters(object)$mode == 'neg', 
                                          'negative mode entries.\n',
                                          'positive mode entries.\n'))
 }
  # n entries
  message('.msp file contains ', prettyNum(sum(grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                                           ignore.case = TRUE)), big.mark = ','), 
          ifelse(Parameters(object)$mode == 'neg', ' negative mode entries.\n',
                 ' positive mode entries.\n'))
  flush.console()
  
  # 2. SMILES calculation
  mfIdxTmp <- mandEntryData$DBSMILES == ''
  
  if(any(mfIdxTmp)){
    message('SMILES code entries were not found in the .msp file.\n')
    flush.console()
    # inchi present and correct mode
    mfIdxTmp <- mfIdxTmp & mandEntryData$DBINCHI != '' 
    mfIdxTmp <- mfIdxTmp & grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                                 ignore.case = TRUE)
    if(any(mfIdxTmp)){
      message(prettyNum(sum(mfIdxTmp), big.mark=','), ' InChI code entries will be converted to SMILES using the command-line interface of OpenBabel (obabel). Please make sure to install this program (openbabel.org) and place in the path...Please wait (processing time highly dependent on number of entries).\n')
      flush.console()
      outDir <- tempfile(pattern = "compMS2Miner")
      dir.create(outDir)
      # create file for openBabel
      in2conv <- mandEntryData$DBINCHI[mfIdxTmp]
      # while loop to determine if every structure converted
      nStrs <- length(in2conv)
      st2conv <- 1
      pb <- txtProgressBar(max=nStrs, style = 3)
      while(st2conv < nStrs){
        inchiFile <- paste0(outDir, '/tmp.inchi')
        smiFile <- paste0(outDir, '/tmp.smi')
        writeLines(paste0(in2conv[st2conv:nStrs], collapse='\n'), con=inchiFile)
        command <- paste0('obabel ', inchiFile, ' -O ', smiFile)
        log <- try(system(command, intern=TRUE, show.output.on.console = FALSE, ignore.stderr = TRUE))
        if(class(log) != "try-error"){
          smiConvTmp <- readLines(smiFile)
          if(length(smiConvTmp) > 0){
            lSmi <- length(smiConvTmp)
            meIdx <- st2conv:{{st2conv - 1} + lSmi}
            mandEntryData$DBSMILES[which(mfIdxTmp)[meIdx]] <- smiConvTmp
            st2conv <- st2conv + lSmi 
          } else {
            #iterate if error
            st2conv <- st2conv + 1  
          }
        }
        setTxtProgressBar(pb, st2conv)
        if(st2conv == nStrs){
          break
        }
      }
      # remove tab delimiter
      mandEntryData$DBSMILES <- gsub('\t', '', mandEntryData$DBSMILES)
      
    }
  } 
  # else if smiles still missing try the inchiKey
  mfIdxTmp <- mandEntryData$DBSMILES == ''
  
  if(all(mfIdxTmp)){
    message('SMILES code entries were still not found in the .msp file. Looking in Comments for InChI codes. Please wait...\n')
    flush.console()
    
    entInfoTmp <- entryInfo[grep('_COMMENTS$', names(entryInfo), ignore.case = TRUE)]
    entryNosTmp <- gsub('_COMMENTS$', '', names(entInfoTmp), ignore.case=TRUE)
    # remove any double entries
    dupIdx <- duplicated(entryNosTmp) == FALSE
    entInfoTmp <- entInfoTmp[dupIdx]
    entryNosTmp <- entryNosTmp[dupIdx]
    indxTmp <- row.names(mandEntryData) %in% entryNosTmp
    # strsplit comments
    inChIsTmp <- sapply(strsplit(entInfoTmp, '"'), function(x){ 
      gsub('.+InChI=', 'InChI=', x[grep('^InChI |^InChI', x, ignore.case=TRUE)], ignore.case = TRUE)})
    mandEntryData$DBINCHI[indxTmp] <- inChIsTmp
    # inchi present and correct mode
    mfIdxTmp <- mfIdxTmp & mandEntryData$DBINCHI != '' 
    mfIdxTmp <- mfIdxTmp & grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                                 ignore.case = TRUE)
    if(any(mfIdxTmp)){
      message(prettyNum(sum(mfIdxTmp), big.mark=','), ' InChI code entries will be converted to SMILES using the command-line interface of OpenBabel (obabel). Please make sure to install this program (openbabel.org) and place in the path...Please wait (processing time highly dependent on number of entries).\n')
      flush.console()
      outDir <- tempfile(pattern = "compMS2Miner")
      dir.create(outDir)
      # create file for openBabel
      in2conv <- mandEntryData$DBINCHI[mfIdxTmp]
      # while loop to determine if every structure converted
      nStrs <- length(in2conv)
      st2conv <- 1
      pb <- txtProgressBar(max=nStrs, style = 3)
      while(st2conv < nStrs){
        inchiFile <- paste0(outDir, '/tmp.inchi')
        smiFile <- paste0(outDir, '/tmp.smi')
        writeLines(paste0(in2conv[st2conv:nStrs], collapse='\n'), con=inchiFile)
        command <- paste0('obabel ', inchiFile, ' -O ', smiFile)
        log <- try(system(command, intern=TRUE, show.output.on.console = FALSE, ignore.stderr = TRUE))
        if(class(log) != "try-error"){
          smiConvTmp <- readLines(smiFile)
          if(length(smiConvTmp) > 0){
          lSmi <- length(smiConvTmp)
          meIdx <- st2conv:{{st2conv - 1} + lSmi}
          mandEntryData$DBSMILES[which(mfIdxTmp)[meIdx]] <- smiConvTmp
          st2conv <- st2conv + lSmi 
          } else {
            #iterate if error
          st2conv <- st2conv + 1  
          }
        }
        setTxtProgressBar(pb, st2conv)
        if(st2conv == nStrs){
          break
        }
      }
      # remove tab delimiter
      mandEntryData$DBSMILES <- gsub('\t', '', mandEntryData$DBSMILES)
      
    }
  } 
  
  if(all(mandEntryData$DBSMILES == '')){
    stop('.msp file does not contain any chemical structure information (SMILES or InChI), therefore no downstream analysis can be performed. Add this information to the entries or try another .msp file.\n')
  }
  
  
  # 3. precursor adduct entry
  mfIdxTmp <- mandEntryData$DBesiType == ''
  # correct mode
  mfIdxTmp <- mfIdxTmp & grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                               ignore.case = TRUE)
  if(any(mfIdxTmp)){
    message('For ', prettyNum(sum(mfIdxTmp), big.mark=','), 
            ' entries no precursor adduct type was identified.\n',
            'The end of each entry name will checked for the adduct type information.\n',
            'In the case that the adduct type is also not detected in the entry name it will be assumed that entries represent the (de-)protonated parent ion depending on polarity. (e.g. [M-H]- or [M+H]+)\n')
    flush.console()
    mandEntryData$DBesiType[mfIdxTmp] <- rep(ifelse(Parameters(object)$mode == 'pos', '[M+H]+', '[M-H]-'), 
                                             sum(mfIdxTmp))
    # see if at end of name
    endNameAdduct <- gsub('.+ ', '', mandEntryData$DBname[mfIdxTmp])
    adductEndIndx <- grepl('M\\+|M\\-', endNameAdduct)
    endNameAdduct[adductEndIndx == FALSE] <- ''
    adductEndIndx[endNameAdduct == mandEntryData$DBname[mfIdxTmp]] <- FALSE
    if(any(adductEndIndx)){
    endNameAdduct <- ifelse(grepl('\\[', endNameAdduct) == FALSE, 
                            paste0('[', endNameAdduct, ']'), endNameAdduct)
    # if necessary add charge 
    chStateIndx <- grepl('\\+$|-$', endNameAdduct)
    if(any(chStateIndx) == FALSE){
      chSym <- ifelse(Parameters(object)$mode == 'neg', '-', '+')
      chState <- paste0(abs(as.numeric(mandEntryData$chState[mfIdxTmp])), chSym) 
      chState <- ifelse(grepl('^0', chState), '', chState)
      chState <- ifelse(chState %in% paste0('1', chSym), '-', chState)
     
      endNameAdduct <- ifelse(chStateIndx == FALSE, paste0(endNameAdduct, chState), endNameAdduct)
      mandEntryData$DBesiType[mfIdxTmp][adductEndIndx] <- endNameAdduct[adductEndIndx]
    }
    } 
  }
  
  # 4. precursor mass calculation if no precursor mass and esi types
  mfIdxTmp <- mandEntryData$DBprecursorMasses == ''
  # correct mode
  mfIdxTmp <- mfIdxTmp & grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                               ignore.case = TRUE)
  if(any(mfIdxTmp)){
    # calculate monoisotopic mass
    message(prettyNum(sum(mfIdxTmp), big.mark=','), ' precursor m/z values were not detected in the .msp file, SMILES codes will therefore be used to calculate the expected monoisotopic mass.\n', 'The ESI adduct type identified will then be used to calculate the expected precursor mass.\n')
    flush.console()
    
    # remove all numbers and punctuation from smiles
    eleComp <- gsub('[[:punct:]]|[[:digit:]]', '', mandEntryData$DBSMILES[mfIdxTmp])
    # capitalize aromatic lower case letters
    for(arom in c('c', 'o', 's', 'n')){
     eleComp <- gsub(arom, toupper(arom), eleComp) 
    }
    # split into elements
    eleComp <- gsub('([[:upper:]])', ' \\1', eleComp)
    # extract mass of most abundant isotope
    # most adundant isotope
    mostAbundIso <-  apply(exactMassEle, 1, function(x){
      monoMassesTmp <- as.numeric(strsplit(x['monoMass'], ' ')[[1]])
      natAbundTmp <- as.numeric(strsplit(x['natAbund'], ' ')[[1]])
      monoMassesTmp <- monoMassesTmp[which.max(natAbundTmp)]})
    names(mostAbundIso) <- exactMassEle$eleSymbol 
    # compute monomass
    monoisotopicMasses <- sapply(strsplit(eleComp, ' '), function(x){
                                 eleCount <- table(x)
                                 mTmp <- mostAbundIso[names(eleCount)] * eleCount
                                 sum(mTmp[!is.na(mTmp)])})
    # catch any 
    esiTypeMassShifts <- data.frame()
    for(add in unique(mandEntryData$DBesiType[mfIdxTmp])){
      mShift <- tryCatch({
        adduct2mass(add)}, error=function(cond){
          data.frame()
        }, warning=function(cond){
          data.frame()})  
      if(nrow(mShift) & all(is.na(mShift[1, ]) == FALSE)){
      esiTypeMassShifts <-  rbind(esiTypeMassShifts, mShift)
      }
    }
    # subset monoisotopic masses
    monoisotopicMasses <- monoisotopicMasses[mandEntryData$DBesiType[mfIdxTmp] %in% esiTypeMassShifts$name]
    # if any were not converted adjust logical index
    mfIdxTmp <- mfIdxTmp & {mandEntryData$DBesiType %in% esiTypeMassShifts$name}
    namesTmp <- esiTypeMassShifts$name
    esiTypeMassShifts <- esiTypeMassShifts$massDiff
    names(esiTypeMassShifts) <- namesTmp
    
    # calc expected mass
    esiTypeMassShifts <- as.numeric(esiTypeMassShifts[mandEntryData$DBesiType[mfIdxTmp]])
    mandEntryData$DBprecursorMasses[mfIdxTmp] <- monoisotopicMasses + esiTypeMassShifts
  } 
  
  # if all wrong mode 
  if(all(grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
               ignore.case = TRUE) == FALSE)){
    stop('None of the .msp file entries were collected in ', 
         ifelse(Parameters(object)$mode == 'neg', 'negative mode.', 'positive mode.'))
  }
  # remove entries where the mandatory entries were not detected
  corrMode <- grepl(Parameters(object)$mode, mandEntryData$DBIonMode, 
                       ignore.case = TRUE)
  missingVals <- apply(mandEntryData[, c('DBSMILES', 'DBprecursorMasses', 'DBesiType')], 1, function(x) all(x != ''))
  # any NAs introduced
  # convert numeric precursor masses
  mandEntryData$DBprecursorMasses <- suppressWarnings(as.numeric(mandEntryData$DBprecursorMasses))
  naIdx <- !is.na(mandEntryData$DBprecursorMasses)
  missingVals <- missingVals & naIdx
  
  if(any(missingVals[corrMode] == FALSE)){
    message(sum(missingVals[corrMode] == FALSE), ' of the ', prettyNum(length(nEntries), big.mark = ','), ' .msp file entries are missing mandatory entry information:\n\n',
            '1. Chemical structural information (either SMILES or InChI).\n',
            '2. A precursor mass (recorded or calculated by compMS2Miner from structure).\n',
            '3. The precursor adduct type.\nand will not be considered.\n')
    flush.console()
  }
    keepEntries <- row.names(mandEntryData)[corrMode & missingVals]
    entryInfo <- entryInfo[gsub('_.+', '', names(entryInfo)) %in% keepEntries]
    massesTmp <- massesTmp[gsub('_.+', '', names(massesTmp)) %in% keepEntries]
    relIntTmp <- relIntTmp[gsub('_.+', '', names(relIntTmp)) %in% keepEntries]
    mandEntryData <- mandEntryData[corrMode & missingVals, , drop=FALSE]
    nEntries <- nEntries[corrMode & missingVals] 
  
  message(prettyNum(length(nEntries), big.mark = ','), ' entries remaining.\n')
  flush.console()
   
  # add database data to table
  DBdf <- cbind(mass=massesTmp, Rel_Intensity=relIntTmp)
  # replace common esi types
  mandEntryData$DBesiType <- ifelse(mandEntryData$DBesiType == '[M+Hac-H]-', "[M-H+CH3COOH]-", mandEntryData$DBesiType)
  mandEntryData$DBesiType <- ifelse(mandEntryData$DBesiType == '[M+FA-H]-', "'[M-H+HCOOH]-'", mandEntryData$DBesiType)
  
  dbNamesTmp <- paste0('DB_', row.names(DBdf))
  message('\ncalculating spectral similarities (dot product >= ', round(minDBDotProdThresh, 2),
          # ' or minimum proportion of spectrum explained >= ', minPropEx, 
          ') between database and composite spectra...\n')
  flush.console()
  # all constituent spectra
  if(verbose == TRUE){ pb <- txtProgressBar(max=length(compSpectra(object)), style = 3)}
  noMatchesV <- as.numeric()
  for(i in 1:length(compSpectra(object))){
    if(verbose==TRUE){setTxtProgressBar(pb, i)}
    massIntTmp <- compSpectra(object)[[i]]
    if(!'Rel_Intensity' %in% colnames(massIntTmp)){
      massIntTmp[, 2] <- {massIntTmp[, 2]/max(massIntTmp[, 2])} * 100
      colnames(massIntTmp)[2] <- 'Rel_Intensity'
    }
    massIntTmp <- massIntTmp[, c('mass', 'Rel_Intensity')]
    ms1MzIndx <- grep('MS1_mz', names(metaData(object)[[i]]), ignore.case=TRUE)
    ms1MzTmp <- unlist(metaData(object)[[i]][ms1MzIndx])[1]
    # calc. ppm diff
    ppmDiff <- ((ms1MzTmp - mandEntryData$DBprecursorMasses)/ ms1MzTmp) * 1E06
    dbMatchMasses <- which(abs(ppmDiff) <= ppmMS1)
    if(length(dbMatchMasses) == 0){
    noMatchesV <- c(noMatchesV, i)  
    next
    } 
    dbMatchEntryIndx <- row.names(mandEntryData)[dbMatchMasses] 
    indxTmp <- row.names(DBdf) %in% dbMatchEntryIndx
    DBdfTmp <- DBdf[indxTmp, , drop=FALSE]
   
    spectrumDBTmp <- rbind(massIntTmp, DBdfTmp)
    specNamesVecTmp <- c(rep('DB_0', nrow(massIntTmp)), dbNamesTmp[indxTmp])
    # padded integer labels
    maxMass <- floor(max(spectrumDBTmp[, 1])) + 10
    labelsTmp <- paste0('(', seq(binSizeMS2, (maxMass - binSizeMS2), binSizeMS2), ',', seq((2 * binSizeMS2), maxMass, binSizeMS2), ']')
    massBinsIndivTmp <- cut(spectrumDBTmp[, 1], breaks=seq(binSizeMS2, maxMass,
                                                           binSizeMS2), 
                            labels=labelsTmp)   
    # empty bins
    indivSpecVec <- tapply(spectrumDBTmp[, 2], paste0(specNamesVecTmp, massBinsIndivTmp), sum)
    # any NAs in names
    indivSpecVec <- indivSpecVec[grepl('NA$', names(indivSpecVec)) == FALSE]
    # identify any absent bins
    allBinNames <- paste0(rep(unique(specNamesVecTmp), each=length(labelsTmp)), rep(labelsTmp, length(unique(specNamesVecTmp))))
    
    # add absent bins as zeros
    allBinsTmp <- rep(0, length(allBinNames))
    names(allBinsTmp) <- allBinNames
    # ensure indivSpecVec is in right order
    allBinsTmp[match(names(indivSpecVec), allBinNames)] <- indivSpecVec
    
    indivSpecMat <- matrix(allBinsTmp, byrow=FALSE, nrow=length(labelsTmp))
    # proportion explained
    propTicEx <- apply(indivSpecMat, 2, function(x){ 
      tmpIndx <- x > 0 
      compSpecIndxTmp <- indivSpecMat[, 1] > 0
      return(sum(indivSpecMat[tmpIndx & compSpecIndxTmp, 1])/ sum(indivSpecMat[compSpecIndxTmp, 1]))})
    # mean all pairwise dotproducts
    # dotProdMat <- t(indivSpecMat) %*% indivSpecMat
    dotProdMat <- crossprod(indivSpecMat)
    sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), nrow=nrow(dotProdMat), 
                            ncol=ncol(dotProdMat), byrow = TRUE) 
    
    dotProdsTmp <- dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
    # first column is the spectrum dot prods
    matchIndxTmp <- dotProdsTmp[, 1] >= minDBDotProdThresh #| propTicEx >= minPropEx
    dotProdsTmp <- dotProdsTmp[matchIndxTmp, , drop=FALSE]
    propTicEx <- propTicEx[matchIndxTmp]
    dotProdsTmp <- dotProdsTmp[-1, 1]
    propTicEx <- propTicEx[-1]
    
    if(length(dotProdsTmp) > 0){
      DBmatchesTmp <- dbMatchEntryIndx[matchIndxTmp[-1]]
      dupIndx <- duplicated(DBmatchesTmp) == FALSE
      DBmatchesTmp <- DBmatchesTmp[dupIndx]
      dotProdsTmp <- dotProdsTmp[dupIndx]
      DBindx <- row.names(mandEntryData) %in% DBmatchesTmp
      compound_msp <- mandEntryData$DBname[DBindx]
      names(compound_msp) <- row.names(mandEntryData)[DBindx]
      names(dotProdsTmp) <- names(compound_msp)
      # if duplicate names then take the highest dot product
      maxDotProd <- tapply(dotProdsTmp, compound_msp, function(x) names(x)[which.max(x)])
      maxDotProd <- names(compound_msp) %in% maxDotProd
      # subset
      DBmatchesTmp <- DBmatchesTmp[maxDotProd]
      dotProdsTmp <- dotProdsTmp[maxDotProd]
      compound_msp <- compound_msp[maxDotProd]
      propTicEx <- propTicEx[maxDotProd]
      DBindx <- row.names(mandEntryData) %in% DBmatchesTmp
      indxTmp <- rownames(DBdf) %in% DBmatchesTmp
      dbMassesInt <- DBdf[indxTmp, , drop=FALSE]
      indxTmp <- gsub('_.+', '', names(entryInfo)) %in% DBmatchesTmp
      entryInfoTmp <- split(entryInfo[indxTmp], gsub('_.+', '', names(entryInfo[indxTmp])))
      
      # if contains smiles then add to the best Annotations
      # if(length(DBSMILES) > 0){
      if(length(DBanno(object)) == 0){
        DBanno(object) <- lapply(1:length(compSpectra(object)), function(x){
          tmpDf <- data.frame(matrix(ncol=10), stringsAsFactors = FALSE)
          colnames(tmpDf) <- c('WebAddress', 'DBid', 'DBname', 'SMILES',
                               'ESI_type', 'SubStr_type', 'observedMass', 
                               'candidateMass', 'ppmMatch', 'BestAnno')
          tmpDf <- tmpDf[-1, , drop=FALSE]
        })
        names(DBanno(object)) <- names(compSpectra(object))
      }
      bestAnnoPrev <- DBanno(object)[[i]]
      bestAnnoSpecDbTmp <- data.frame(matrix('', ncol=ncol(bestAnnoPrev), 
                                             nrow=length(DBmatchesTmp)), 
                                      stringsAsFactors = FALSE) 
      colnames(bestAnnoSpecDbTmp) <- colnames(bestAnnoPrev)
      bestAnnoSpecDbTmp$WebAddress <- rep(basename(mspFile), length(DBmatchesTmp))
      bestAnnoSpecDbTmp$DBid <- DBmatchesTmp
      bestAnnoSpecDbTmp$DBname <- gsub(';.+', '', compound_msp)
      bestAnnoSpecDbTmp$SMILES <- gsub('\\*', '', mandEntryData$DBSMILES[DBindx])
      bestAnnoSpecDbTmp$ESI_type <- mandEntryData$DBesiType[DBindx]
      bestAnnoSpecDbTmp$SubStr_type <- rep('', nrow(bestAnnoSpecDbTmp))
      bestAnnoSpecDbTmp$observedMass <- rep(ms1MzTmp, nrow(bestAnnoSpecDbTmp))
      bestAnnoSpecDbTmp$candidateMass <- mandEntryData$DBprecursorMasses[DBindx]
      bestAnnoSpecDbTmp$ppmMatch <- ppmDiff[DBindx]
      bestAnnoSpecDbTmp$BestAnno <- rep(TRUE, nrow(bestAnnoSpecDbTmp))
      bestAnnoSpecDbTmp <- bestAnnoSpecDbTmp[grepl('not_identified', bestAnnoSpecDbTmp$DBname) == FALSE, , drop=FALSE]
      if(nrow(bestAnnoSpecDbTmp) > 0){
      # rbind
      bestAnnoPrev <- rbind(bestAnnoSpecDbTmp, bestAnnoPrev) 
      bestAnnoPrev <- bestAnnoPrev[duplicated(paste0(bestAnnoPrev$WebAddress, bestAnnoPrev$DBid)) == FALSE, , drop=FALSE]
      DBanno(object)[[i]] <- bestAnnoPrev
      }
      # add assignments to comments table 
      if(nrow(bestAnnoSpecDbTmp) > 0){
      if(autoPossId ==  TRUE){
            # met Id comments table
            metIDcomments <- Comments(object)
            # add possible annotations to metIDcomments
            dpIndx <- which.max(dotProdsTmp)
            indxTmp <- metIDcomments$compSpectrum %in% names(compSpectra(object))[i]
            metIDcomments$possible_identity[indxTmp] <- bestAnnoSpecDbTmp$DBname[dpIndx] 
            metIDcomments$ESI_type[indxTmp] <- bestAnnoSpecDbTmp$ESI_type[dpIndx]
            metIDcomments$user_comments[indxTmp] <- 'metID.matchSpectralDB'
            Comments(object) <- metIDcomments
          } # end autoPossId
      }
      # }
      compound_msp <- paste0(compound_msp, '__',
                             gsub('\\.msp$', '', basename(mspFile)))
      names(entryInfoTmp) <- compound_msp
      compound_msp <- unlist(mapply(rep, compound_msp, each=table(row.names(dbMassesInt)), SIMPLIFY = FALSE))
      
      dotProdsTmp <- unlist(mapply(rep, dotProdsTmp, each=table(row.names(dbMassesInt)), SIMPLIFY = FALSE))
      propTicEx <- unlist(mapply(rep, propTicEx, each=table(row.names(dbMassesInt)), SIMPLIFY = FALSE))
      dbMassesInt <- cbind(dbMassesInt, compound_msp, dotProductScore=dotProdsTmp, propTIC_explained=propTicEx)
      if(length(spectralDB(object)[[i]]) > 0){
      spectralDB(object)[[i]]$dbSpectra <- rbind(spectralDB(object)[[i]]$dbSpectra, 
                                                dbMassesInt)
      spectralDB(object)[[i]]$entryInfo <- c(spectralDB(object)[[i]]$entryInfo, 
                                            entryInfoTmp)
      } else {
      spectralDB(object)[[i]] <- list(dbSpectra=dbMassesInt, entryInfo=entryInfoTmp)
      }
      # message(i, '\n', max(dotProdsTmp[2:nrow(dotProdsTmp), 1]), '\n')
    } else {
      noMatchesV <- c(noMatchesV, i)
    }
  }
  # set the difference between those with any previous matches and new no matches vector
  noMatchesV <- setdiff(noMatchesV, which(sapply(spectralDB(object), length) > 0))
  
  nMatched <- length(compSpectra(object))-length(noMatchesV)
  message('\n', nMatched, ' composite spectra (', 
          round(((length(compSpectra(object))-length(noMatchesV))/length(compSpectra(object))) * 100, 1), 
          '%) currently matched to spectral databases.\n')
  flush.console()
  if(nMatched > 0){
  message('These were added to the Comments table and flagged as metID.matchSpectralDB.\n')
  flush.console()
  }
  return(object)
}) # end function
