#' cfm fragment graph generation from table of annotations
#'@param bestAnnoSubRow unique compMS2@BestAnno entries (only M-H (neg mode) and M+H (pos mode) can be in silico fragmented by CFM and no SubStr_types).
#'@param fragGraphGenExe character full path to fraggraph-gen.exe file (internal to compMS2Miner package).
#'@param compSpecAll data.frame 3 columns mass, intensity and comp spectrum index number. 
#'@param mode character ionization polarity (either 'pos' or 'neg').
#'@param frag_mzabs numeric delta predicted-observed fragment mass accuracy for matching.
#'@return if fraggraph-gen process completed then a list of fragments matched to corresponding composite spectra are return
cfmFragGraphGen <- function(bestAnnoSubRow=NULL, fragGraphGenExe=NULL, compSpecAll=NULL, keepTempFiles=FALSE, mode='pos', frag_mzabs=0.05){
  if(keepTempFiles == TRUE){
    tmpCFMdir <-  paste0(getwd(), '/CFM_results/')
    suppressWarnings(dir.create(tmpCFMdir))
    dir.create(paste0(tmpCFMdir, bestAnnoSubRow$DBid))
    csvNameTmp <- paste0(tmpCFMdir, bestAnnoSubRow$DBid, "/cfmFragOutput.csv")  
    fragMode <- ifelse(mode == 'pos', ' 1 + ', ' 1 - ')
  } else {
    tmpCFMdir <-  tempfile(pattern = "compMS2Miner")
    dir.create(tmpCFMdir)
    fragMode <- ifelse(mode == 'pos', ' 1 + ', ' 1 - ')
    csvNameTmp <- paste0(tmpCFMdir, "/tempCFMoutput", bestAnnoSubRow$DBid, ".csv")   
  }
  commandTmp <- paste0('"', fragGraphGenExe, '" "',  bestAnnoSubRow$SMILES, '"', fragMode, ' fragonly "', csvNameTmp,'"')  
  cmdRes <- system(commandTmp, intern = TRUE)
  
  nonZeroSize <- file.info(csvNameTmp)$size > 0
  # load csv if necc
  if(nonZeroSize){
    resTmp <- readLines(csvNameTmp)
    resTmp <- do.call(rbind, strsplit(resTmp, ' '))
    specIdsTmp <- as.numeric(strsplit(as.character(bestAnnoSubRow$querySpecDbId), ' ')[[1]])
    compSpecSubTmp <- compSpecAll[compSpecAll$featureSubSet %in% specIdsTmp, , drop=FALSE]
    sumTic <- tapply(compSpecSubTmp$intensity, compSpecSubTmp$featureSubSet, sum)  
    indxTmp <- match(compSpecSubTmp$featureSubSet, names(sumTic))
    compSpecSubTmp$sumTIC <- sumTic[indxTmp]
    matchedPeaks <- as.numeric()
    for(j in 1:nrow(resTmp)){
      # calc mass difference
      indxTmp <- which(abs(compSpecSubTmp$mass - as.numeric(resTmp[j, 2])) < frag_mzabs)
      if(length(indxTmp) > 0){
        names(indxTmp) <- rep(j, length(indxTmp))
        matchedPeaks <- c(matchedPeaks, indxTmp)
      }
    }
    
    if(length(matchedPeaks) > 0){
      compSpecSubTmp <- compSpecSubTmp[matchedPeaks, , drop=FALSE]
      compSpecSubTmp <- cbind(compSpecSubTmp, resTmp[as.numeric(names(matchedPeaks)), , drop=FALSE])
      featureSubSetTmp <- compSpecSubTmp$featureSubSet
      compSpecSubTmp$featureSubSet <- NULL
      colnames(compSpecSubTmp)[4:ncol(compSpecSubTmp)] <- c('CFM_rank', 'CFM_mass', 'CFM_fragSMILES')
      compSpecSubTmp$CFM_rank <- paste0(compSpecSubTmp$CFM_rank, ' of ', nrow(resTmp), ' fragments')
      compSpecSubTmp$CFM_fragPropEx <- compSpecSubTmp$intensity/compSpecSubTmp$sumTIC
      maxFragIntCompSpec <- paste0(featureSubSetTmp, '_', compSpecSubTmp$mass)
      maxFragIntCompSpec <- tapply(compSpecSubTmp$CFM_fragPropEx, maxFragIntCompSpec, max)
      sumExTmp <- tapply(maxFragIntCompSpec, gsub('_.+', '', names(maxFragIntCompSpec)), sum)
      
      indxTmp <- match(featureSubSetTmp, names(sumExTmp))
      compSpecSubTmp$CFM_totPropEx <- sumExTmp[indxTmp]
      # add in best Anno info
      compSpecSubTmp <- cbind(compSpecSubTmp, bestAnnoSubRow[rep(1, nrow(compSpecSubTmp)), c('WebAddress', 'DBid', 'DBname', 'SMILES'), drop=FALSE])
      compSpecSubTmp <- split(compSpecSubTmp, f=featureSubSetTmp)
      return(compSpecSubTmp)
    } # if any matches
  } # if non-zero file
} # end cfmFragGraphGen function
