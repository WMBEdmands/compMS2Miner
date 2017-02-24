#' automatically annotate compound classes comments table
#' @param object a compMS2 class object.
#' @return the compMS2 class object with compound classes (currently only phospholipids supported) contained in
#' data(lipidAbbrev) annotated and interpreted lipid tail information added to the Comments table
#' accessible by Comments(object).
#' @param overWrite logical should any existing compound_class types be overwritten in the
#' comments table (default = FALSE). 
#' @param minSimScore numeric mimimum chemical similarity score (default = 0.8),
#' all chemical fingerprints with similarity above this score to entries of HMDB,
#' drugBank and T3DB will be considered for the compound class. The highest 
#' chemical similarity score above this threshold will be considered as the 
#' predicted compounds class. In the event of multiple tied maximum chemical
#' similarity scores the most frequent compound class will added to the comments
#' table.
#' @details This function attempts to automatically add chemical taxonomy information
#' to the metID comments table accessible using Comments(object). Pre-existing
#' compound class information can be overwritten. Additionally a predicted
#' compound class is added for any compound not contained in HMDB, drugBank or
#' T3DB above a minimum chemical similarity score.
#' @return a compMS2 class object with the compound class of the identified
#' compound added to the metID comments table. Comments(object). Additionally
#' a barplot is generated summarizing the identified compound class information
#' returned.  
#' @export
setGeneric("annoCompoundClass", function(object, ...) standardGeneric("annoCompoundClass"))

setMethod("annoCompoundClass", signature = "compMS2", function(object, overWrite=FALSE, minSimScore=0.8, ...){
  require(ChemmineR)
metIDcomments <- Comments(object)  
if(overWrite ==  TRUE){
  metIDcomments$compound_class <- ''
}
lysoIndx <- grepl('/[0]:[0]\\)|\\([0]:[0]/', metIDcomments$possible_identity)
for(i in 1:nrow(lipidAbbrev)){
  indxTmp <- grep(lipidAbbrev$regexpr[i], metIDcomments$possible_identity, 
                  ignore.case = TRUE)
  if(length(indxTmp) > 0){
    metIDcomments$compound_class[indxTmp] <- gsub('.+ ', '', lipidAbbrev$Abbreviation[i])
    # add fatty acid tail type
    fattyAcidType <- gsub('[A-z]|\\-', '', metIDcomments$possible_identity[indxTmp])
    # split based on parentheses
    fattyAcidType <- strsplit(fattyAcidType, '\\(|\\)|/')
    # generate tail info
    fattyAcidType <- sapply(fattyAcidType, function(x){
      x <- x[x != '']
      nTails <- grep('\\:', x)
      saturation <- gsub('.+\\:', '', x[nTails])
      nCarb <- gsub('\\:.+', '', x[nTails])
      tmpIndx <- {nCarb == '0' & saturation == '0'} == FALSE
      saturation <- saturation[tmpIndx]
      nCarb <- nCarb[tmpIndx]
      nTails <- nTails[tmpIndx]
      fAType <- vector('character', length(nTails))
      for(j in 1:length(fAType)){
      if(saturation[j] == '0'){
        fAType[j] <- 'saturated'
      } else if(saturation[j] == '1'){
        fAType[j] <- 'MUFA'
      } else {
        fAType[j] <- 'PUFA'
      }
      }
      fAType <- paste0('fatty acid tail(s): ' ,
                       paste0('(', 1:length(nTails), '). ', nCarb, ' carbon chain - ', fAType, collapse=' '))
    })
    # add info to user comments
    # replace any existing comments
    metIDcomments$user_comments[indxTmp] <- gsub('; fatty acid tail\\(s\\):.+|fatty acid tail\\(s\\):.+', '', metIDcomments$user_comments[indxTmp]) 
    metIDcomments$user_comments[indxTmp] <- paste0(metIDcomments$user_comments[indxTmp], ifelse(metIDcomments$user_comments[indxTmp] == '', '', '; '),
                                                   fattyAcidType)
  }
}
metIDcomments$compound_class[lysoIndx] <- paste0('lyso', metIDcomments$compound_class[lysoIndx])

# compound class based on HMDB, drugBank and T3DB
# index which bestAnnos have entries
emptyIndx <- sapply(BestAnno(object), nrow) > 0
if(all(emptyIndx == FALSE)){
  stop('none of composite spectra have any dbAnnotate matches')
}

smilesDf <- do.call(rbind, BestAnno(object)[emptyIndx])
nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
smilesDf$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[emptyIndx]),
                                           each=nrowBestAnno))

# subset
smilesDf <- smilesDf[, colnames(smilesDf) %in% c('DBid', 'SMILES', 'DBname',
                                                 'SubStr_type', 'ESI_type', 'WebAddress',
                                                 'chemFP', 'specNamesTmp')]
# match comments to smiles data.frame
idxTmp <- match(paste0(metIDcomments$compSpectrum, '_', metIDcomments$possible_identity),
                paste0(smilesDf$specNamesTmp, '_', smilesDf$DBname))
# if already compound class then do not consider
alreadyCompClass <- 1:nrow(metIDcomments) %in% which(!is.na(idxTmp)) & metIDcomments$compound_class != ''
idxTmp[alreadyCompClass] <- NA
smilesDf <- smilesDf[idxTmp[!is.na(idxTmp)], , drop=FALSE]
if(('chemFP' %in% colnames(smilesDf)) == FALSE){
  # add new column
  smilesDf$chemFP <- ''
  # if possible and bitsChemFP equal to 1024 match as many ids as possible to internal
  # databases
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
    indxTmp <- match(smilesDf$DBid, allDbIds)
  smilesDf$chemFP[!is.na(indxTmp)] <- names(allDbIds)[indxTmp[!is.na(indxTmp)]]
  # empty chemFP
  indxTmp <- smilesDf$chemFP == ''
  if(any(indxTmp)){  
    # convert first to SMILES set object and then SDF (ChemmineR/OB)
    message('Converting ', sum(indxTmp), ' SMILES codes -> SDF -> atom-pair descriptors -> 1024 most common atom pairs (DrugBank.ca) chemical fingerprints this only needs to be performed once...please wait\n')
    flush.console()
    SDFtmp <- tryCatch({
      suppressWarnings(ChemmineR::smiles2sdf(gsub('\\*', '', smilesDf$SMILES[indxTmp])))
    }, error=function(cond) {
      return(0)
    }, warning=function(cond){
      return(0)})
    if(!is.numeric(SDFtmp)){
      SDFtmp@ID <- smilesDf$DBid[indxTmp]
      
      # if necessary alert user to failed apset converts
      errorIndx <- ChemmineR::validSDF(SDFtmp)
      if(any(errorIndx == FALSE)){
        message('The following database entries could not be converted to atom-pair ',
                'descriptors:\n', paste0(smilesDf$DBid[indxTmp][errorIndx == FALSE], ' ', paste0('(', smilesDf$DBname[indxTmp][errorIndx == FALSE], ')'), sep='\n'),
                'Chemical similarity scores cannot be calculated.')
        flush.console()
        # remove invalid sdf
        SDFtmp <- SDFtmp[errorIndx]
        # remove from smilesDf
        smilesDf <- smilesDf[-which(smilesDf$DBid %in% smilesDf$DBid[indxTmp][errorIndx == FALSE]), , drop=FALSE]
        smilesDf <- smilesDf[-which(indxTmp)[errorIndx == FALSE], , drop=FALSE]
        indxTmp <- smilesDf$chemFP == ''
      }
      # create APset atom pair descriptors files
      #   tryCatch({
      apsetTmp <- suppressWarnings(ChemmineR::sdf2ap(SDFtmp))
      
      # create n bit finger prints of apset using ChemmineR 
      if(length(apsetTmp) > 0){
        fpsetTmp  <-  ChemmineR::desc2fp(apsetTmp,  descnames=1024,
                                         type="matrix")
        
        smilesDf$chemFP[indxTmp] <- apply(fpsetTmp, 1, function(x) 
          paste0(which(x == 1), collapse=';'))
      }
    }
  }
}
  # see if database ids are present
  availDbs <- c('HMDB', 'drugBank', 'T3DB')
    # attach 
    allDbIds <- do.call(c, lapply(availDbs, function(x) eval(parse(text=paste0('data(', x, ')')))))
    # extract dbids
    allDbIds <- do.call(c, lapply(availDbs, function(x) eval(parse(text=paste0(x, '$Unique_DB_ID')))))
    # extract chemical fingerprint
    names(allDbIds) <- do.call(c, lapply(availDbs, function(x){
      eval(parse(text=paste0(x, '$compound_class')))}))
    # extract chemical fingerprint
    allDbChemFP <- do.call(c, lapply(availDbs, function(x){
      eval(parse(text=paste0(x, '$chemFP')))}))
    # convert matrix
    allDbChemFP <-  do.call(rbind, lapply(allDbChemFP, function(x){
      indxTmp <- as.numeric(strsplit(x, ';')[[1]])
      return(ifelse(1:1024 %in% indxTmp, 1, 0))}))
    row.names(allDbChemFP) <- allDbIds
    # if conjugates add these in from db id
    smilesDf$compound_class <- ifelse(grepl('.+_', smilesDf$DBid), 
                                      paste0(' (', gsub('.+_', '', smilesDf$DBid), ' conjugate)'), '')
  classKnownIdx <- match(gsub('_.+', '', smilesDf$DBid), allDbIds)
  if(any(!is.na(classKnownIdx))){
    smilesDf$compound_class[!is.na(classKnownIdx)] <- paste0(names(allDbIds)[classKnownIdx[!is.na(classKnownIdx)]], 
                                                             smilesDf$compound_class[!is.na(classKnownIdx)])
  } 
  if(any(is.na(classKnownIdx))){
    cat('calculating chemical similarity between', sum(is.na(classKnownIdx)), 'compounds and entries of HMDB, DrugBank and T3DB. The compound class for these compounds in the metID comments table will be predicted.\n')
    simIdx <- which(is.na(classKnownIdx))
    for(j in simIdx){
    quMol <- vector('numeric', 1024)
    quMol[as.numeric(strsplit(smilesDf$chemFP[j], ';')[[1]])] <- 1
    simMol <- ChemmineR::fpSim(quMol, allDbChemFP, sorted = FALSE)
    maxScoreIdx <- round(simMol, 2) %in% round(max(simMol), 2)
    aboveScoreIdx <- simMol >= minSimScore
    if(any(aboveScoreIdx)){
    maxSimTmp <- table(names(allDbIds)[maxScoreIdx])
    smilesDf$compound_class[j] <- paste0(names(maxSimTmp)[which.max(maxSimTmp)], 
                                         smilesDf$compound_class[j], 
                                         ' (predicted by annoCompoundClass)')  
    }
    }
  }
  metIDcomments$compound_class[!is.na(idxTmp)] <- smilesDf$compound_class
  # add contaminants and no annotations
  contamIdx <- apply(metIDcomments, 1, function(x) any(grepl('contaminant', x, ignore.case = TRUE)))
  metIDcomments$compound_class[contamIdx] <- 'possible contaminant'
  noAnnoIdx <- metIDcomments$possible_identity == 'no_annotations'
  metIDcomments$compound_class[noAnnoIdx] <- 'no_annotations'
  notKnownIdx <-  metIDcomments$possible_identity == '' | metIDcomments$possible_identity == 'unclear'
  metIDcomments$compound_class[notKnownIdx] <- 'unknown'
  # make compound class barplot
  predIdx <- grepl(' \\(predicted by annoCompoundClass\\)', metIDcomments$compound_class)
  compoundClasses <- gsub(' \\(predicted by annoCompoundClass\\)', '', metIDcomments$compound_class)
  # rename for consistency
  namesTmp <- sapply(compoundClasses, function(x) compoundClasses[grep(paste0('^', x, '$'), compoundClasses, ignore.case=TRUE)[1]])
  namesTmp <- ifelse(is.na(namesTmp), names(namesTmp), namesTmp)
  metIDcomments$compound_class <- ifelse(predIdx, paste0(namesTmp, ' (predicted by annoCompoundClass)'), namesTmp)
  namesTmp <- namesTmp[namesTmp != '']
  compClassFreq <- table(namesTmp)
 
  par(mar=c(18,4.1,4.1,2.1))
  ylim <- c(0, 1.1 * max(compClassFreq))
  sumClasses <- length(compClassFreq[grepl('unknown|no_annotations|possible contaminant', names(compClassFreq)) ==  FALSE])
  xx <- barplot(as.numeric(compClassFreq), xaxt = 'n', xlab = '', width = 0.85, ylim = ylim,
                ylab = "Frequency", las=2, col='lightblue', 
                main=paste0(length(compClassFreq), ' compound classes (n=', sum(compClassFreq), ' total)'))
  ## Add text at top of bars
  text(x = xx, y = compClassFreq, label = compClassFreq, pos = 3, cex = 0.8, col = "red")
  ## Add x-axis labels 
  axis(1, at=xx, labels=substr(names(compClassFreq), 1, 50), tick=FALSE, las=2, 
       line=-0.5, cex.axis=0.9)
# add comments back to object
Comments(object) <- metIDcomments
return(object)
}) # end function
