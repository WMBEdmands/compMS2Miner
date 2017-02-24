#' metFrag command line function for localSDF files 
#' 
#' @param massSpectrum data.frame composite spectrum consisting of two columns mass and intensity.
#' @param precMass numeric the MS1 m/z (precursor mass).
#' @param compSpecName character name of composite spectrum for directory and file naming.
#' @param dbEntryTable data.frame with the requisite information for the SDFtmp localSDF database consisting of at least 4 columns named 1. 'WebAddress', 2. 'DBid', 3. 'DBname', 4. 'SMILES'.
#' @param metFragJar character full path to metFragCL.jar file (extdata in compMS2Miner package).
#' @param SDFtmp an "SDFset" class object of SDF file for the localSDF database search of metFragCL.
#' @param keepTempFiles logical default = FALSE, sdf, mf and results files will
#' be created as temporary files otherwise if TRUE files will be retained in named subdirectories (see details).
#' @param mode character ionization polarity (either 'pos' or 'neg').
#' @param frag_mzabs numeric delta predicted-observed fragment mass accuracy for matching.
#' @param esiList named numeric vector of electrospray type numbers for metFrag params file. e.g. positive mode 
#' \tabular{llll}{
#' \cr M+H \tab M+NH4  \tab M+Na   \tab M+K 
#' \cr 1    \tab 18    \tab 23    \tab 39
#' }
#' @param maxTreeDepth numeric fragments of fragments? (default = 1 i.e. only direct daughter ions of the structure will be considered). Setting the tree depth to higher values may cause the metFragCL to take longer.
#' @source \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/} developed based on the command line .jar file (MetFrag2.3-CL.jar) downloaded on 2016/07/12
#' \enumerate{
#' \item MetFrag relaunched: incorporating strategies beyond in silico fragmentation: C Ruttkies, E L Schymanski, S Wolf, J Hollender, S Neumann Journal of Cheminformatics 2016 8:3
#' \item In silico fragmentation for computer assisted identification of metabolite mass spectra: S Wolf, S Schmidt, M Müller-Hannemann, S Neumann BMC bioinformatics 11 (1), 148
#' }
#' @return if MetFrag2.3-CL.jar process completed then a data.frame containing any fragments matched to the composite mass spectra are returned.
metFragCl <- function(massSpectrum=NULL, precMass=NULL, compSpecName=NULL,
                      dbEntryTable=NULL, metFragJar=NULL, SDFtmp=NULL, 
                      keepTempFiles=FALSE, mode='pos', frag_mzabs=0.05, esiList=NULL,
                      maxTreeDepth=1){
  # error handling
  stopifnot(is.data.frame(massSpectrum))
  stopifnot(is.numeric(precMass))
  stopifnot(is.character(compSpecName))
  stopifnot(file.exists(metFragJar))
  stopifnot(class(SDFtmp) == 'SDFset')
  
  if(keepTempFiles == TRUE){
    tmpMetFragDir <-  paste0(getwd(), '/MetFrag_results/')
    suppressWarnings(dir.create(tmpMetFragDir))
    tmpMetFragDir <- paste0(tmpMetFragDir, compSpecName, '/')
    suppressWarnings(dir.create(tmpMetFragDir))
    peakListFileName <- paste0(tmpMetFragDir, compSpecName, 'peakList.txt')
    paramFileName <- paste0(tmpMetFragDir, compSpecName, 'metFragParam.txt')
    localSDFNameTmp <- paste0(tmpMetFragDir, compSpecName, 'localSDF.sdf')
    outputFileName <- paste0(tmpMetFragDir, compSpecName, '.csv')
  } else {
    tmpMetFragDir <-  paste0(tempfile(pattern = "compMS2Miner"), '/')
    suppressWarnings(dir.create(tmpMetFragDir))
    peakListFileName <- paste0(tmpMetFragDir, compSpecName, 'peakList.txt')
    paramFileName <- paste0(tmpMetFragDir, compSpecName, 'metFragParam.txt')
    localSDFNameTmp <- paste0(tmpMetFragDir, compSpecName, 'localSDF.sdf')
    outputFileName <- paste0(tmpMetFragDir, compSpecName, '.csv')
  }
  # results 
  resultsTmp <- data.frame()
  peakListFile <- writeLines(paste(paste(massSpectrum$mass, massSpectrum$intensity),
                                   collapse="\n"), con=peakListFileName)
 
  for(esiTypeTmp in unique(dbEntryTable$ESI_type)){
    nMetFragOutRows <- 0
  SDFSubTmp <- SDFtmp[dbEntryTable$ESI_type %in% esiTypeTmp]
  if(length(SDFSubTmp) > 0){
 
  
  while(nMetFragOutRows < length(SDFSubTmp)){
  
  ChemmineR::write.SDF(SDFSubTmp[(nMetFragOutRows + 1):length(SDFSubTmp)], localSDFNameTmp)
  
  paramFileTmp <- paste0('PeakListPath = ', peakListFileName, '\n',
                         'MetFragDatabaseType = LocalSDF', '\n',
                         'LocalDatabasePath = ', localSDFNameTmp, '\n',
                         'IonizedPrecursorMass = ', precMass, '\n',
                         'FragmentPeakMatchAbsoluteMassDeviation = ', frag_mzabs, '\n',
                         'FragmentPeakMatchRelativeMassDeviation = ', 5, '\n',	
                         'PrecursorIonMode = ', esiList[names(esiList) %in% esiTypeTmp], '\n',	
                         'IsPositiveIonMode = ', ifelse(mode == 'pos', 'True', 'False'), '\n',
                         'MetFragScoreTypes = FragmenterScore\n',	
                         'MetFragScoreWeights = 1.0\n',
                         'MetFragCandidateWriter = CSV\n',	
                         'SampleName = ', compSpecName, '\n',	
                         'ResultsPath = ', tmpMetFragDir, '\n',
                         'MaximumTreeDepth = ', maxTreeDepth, '\n')#,
                        #'MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter%\n',	
                        # 'MetFragPostProcessingCandidateFilter = InChIKeyFilter\n')	
  # write param file                       
  writeLines(paramFileTmp, con=paramFileName)
  
  command <- paste0('java -jar "', metFragJar,'" "', paramFileName, '"')
  log <- system(command, intern=TRUE, show.output.on.console = FALSE, ignore.stderr = TRUE)
  
  tmpMfOut <- tryCatch({
   read.csv(paste0(tmpMetFragDir, compSpecName, '.csv'), header=TRUE, stringsAsFactors = FALSE)}, error=function(cond) {
      data.frame()
    }, warning=function(cond){
      data.frame()})
  # add unique db id
  if(nrow(tmpMfOut) > 0){
  tmpMfOut$DBid <- SDFSubTmp@ID[tmpMfOut$Identifier + nMetFragOutRows]
  if(!exists('metFragOutput')){
    metFragOutput <- data.frame(matrix(NA, ncol=ncol(tmpMfOut), nrow = length(SDFSubTmp)), 
                                  stringsAsFactors = FALSE) 
    colnames(metFragOutput) <- colnames(tmpMfOut)
  }
  indxTmp <- {nMetFragOutRows + 1}:{nrow(tmpMfOut) + nMetFragOutRows}
  metFragOutput[indxTmp, ] <- tmpMfOut
  } else if(nrow(tmpMfOut) == 0 & length(SDFSubTmp) == 1){
  metFragOutput <- data.frame(result='?metID.metFrag returned no results', stringsAsFactors = FALSE)
  break
  }
  nMetFragOutRows <- nMetFragOutRows + {nrow(tmpMfOut)  + 1}
  }
  
    if(!exists('metFragOutput')){
      metFragOutput <- data.frame(result='?metID.metFrag returned no results', stringsAsFactors = FALSE)  
    }
  # if not results returned
  if(ncol(metFragOutput) > 1){
  # error indx
  errIndx <- apply(metFragOutput, 1, function(x) all(is.na(x)))
  # add dbId to any errors
  metFragOutput$DBid[errIndx] <- SDFSubTmp@ID[errIndx]
  if(nrow(metFragOutput) > 0){
    # recalculate relative mf score
    metFragOutput$Score <- metFragOutput$FragmenterScore/ max(metFragOutput$FragmenterScore)
    dbEntryTableSub <- dbEntryTable[match(metFragOutput$DBid, dbEntryTable$DBid), c('WebAddress', 'DBname', 'SMILES', 'ESI_type'), drop=FALSE]
    # mass and intensity values
    massIntTmp <- do.call(rbind, lapply(strsplit(as.character(metFragOutput$FormulasOfExplPeaks), ':|;'),
              function(x){
              if(!is.na(x[1])){
              massTmp <-  round(as.numeric(x[seq(1, length(x), 2)]), 4)
              # add precursor
              massTmp <- c(massTmp, massSpectrum$mass[which(abs(massSpectrum$mass - precMass) < frag_mzabs)])
              intTmp <- round(massSpectrum$intensity[round(massSpectrum$mass, 4) %in% round(massTmp, 4)], 2)
              fragSmiTmp <- x[seq(2, length(x), 2)]
              propIntTmp <- round(sum(intTmp)/ sum(massSpectrum$intensity), 3)
              return(c(propIntTmp,
                       paste0(massTmp, collapse = '; '), 
                       paste0(intTmp, collapse = '; '), 
                       paste0(fragSmiTmp, collapse = '; ')))
              } else {
              return(rep(NA, 4))
              }}))
      
    colnames(massIntTmp) <- c('propIntEx', 'mass', 'intensity', 'frag_smiles')
    
    metFragOutput <- metFragOutput[, c('DBid', 'Score', 'FragmenterScore_Values',
                                       'FragmenterScore', 'NoExplPeaks', 
                                       'NumberPeaksUsed'), drop=FALSE]
    metFragOutput <- cbind(dbEntryTableSub, metFragOutput, massIntTmp)
    resultsTmp <- rbind(resultsTmp, metFragOutput)
    rm(metFragOutput)
  }
  } else {
  rm(metFragOutput)  
  }
  } 
  }

  if(nrow(resultsTmp) == 0){
  return(data.frame(result='?metID.metFrag returned no results', stringsAsFactors = FALSE))
  } else {
  return(resultsTmp)
  }
  # if more than zero rows returned
} # end metFragCl function
