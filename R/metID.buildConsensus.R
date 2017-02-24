#' Build consensus combinatorial metabolite identification
#' 
#' @description this function is designed to be used at the end of combinatorial
#' metabolite identification process. It evaluates the multiple layers of evidence
#' which are currently accumulated in the CompMS2 class object to automatically  
#' rank possible annotations and identify the annotation with the greatest weight
#' of evidence for every composite spectrum. 
#' 
#' @param object a "compMS2" class object.
#' @param include character vector of 6 options to build consensus combinatorial
#' metabolite identification see Details below for a description of each. If
#' specific options are not supplied as a character vector then the default 
#' is to consider all 7. i.e. 
#' c('massAccuracy', 'spectralDB', 'inSilico', 'rtPred', 'chemSim', 'pubMed', 
#' 'substructure').
#' @param metIDWeights numeric vector equal in length to include vector (see above).
#' Default is NULL and a simple arithmetic mean will be calculated for all the 
#' metabolite identification options included. The metIDWeights will be used to
#' calculate a weighted mean of the combination of metabolite identification options. 
#' This option can be used to generate a custom metabolite identification setting
#' which best annotates the unknown metabolites. \emph{N.B.} The sum of the 
#' metIDWeights vector must be 1. e.g. include= c('massAccuracy', 'spectralDB',
#' 'inSilico') and metIDWeights=c(0.2, 0.5, 0.3) therefore massAccuracy will be
#' given a weight of 0.2 (20\%), spectralDB matches will be given a weight of 0.5
#' (50\%) and \emph{in silico} fragmentation score will be given a weight of 
#' 0.3 (30\%). rtPred (predicted retention time), chemSim (nearest neighbour chemical
#' similarity score) and pubMed (number of pubmed citations) will not be included. 
#' @param autoPossId logical if TRUE the function will automatically add the names
#'  of the top annotation based on mean consensus annotation score to the 
#'  "metID comments" table (default = FALSE). Caution if TRUE this will overwrite 
#'  any existing possible_identities in the "metID comments" table. 
#'  This functionality is intended as an automatic metabolite annotation 
#'  identification tool prior to thorough examination of the data in 
#'  \code{\link{compMS2Explorer}} as part of an objective and seamless 
#'  first-pass annotation workflow. The mean build consensus score can consist
#'  of many orthogonal measurements of metabolite identification and a means
#'  to rapidly rank metabolite annotations.
#' @param minMeanBCscore numeric minimum mean consensus score (values between 0-1),
#'  if argument autoPossId is TRUE any metabolite annotations above this value
#'  will be automatically added to the "metID comments" table. (if argument not
#'  supplied the default is the upper interquartile range of the mean BC score).
#' @param possContam numeric how many times does a possible annotation have
#' to appear in the automatically generated possible annotations for it to be
#' considered a contaminant and therefore not added to the "metID comment" table (default = 3, i.e. if a database name appears more than 3 times in the 
#'  automatic annotation table it will be removed).
#' @param verbose logical if TRUE display progress bars.
#'@details Specifically the function looks at the following 7 pieces of 
#' evidence:
#' \enumerate{
#' \item \strong{"massAccuracy"} monoisotopic mass similarity. Absolute mass similarity between 0 and 
#' the upper mass accuracy limit (default 10 ppm) are used to generate a ranking
#' score between 0-1.
#' \item \strong{"spectralDB"} spectral database match. If a match has been made to a spectral database
#' using the function \code{\link{metID.matchSpectralDB}} then a combination of
#' the dot product score and proportion of the composite spectrum explained is
#' used to rank the annotations. A score is determined between 0-1 based on the
#' average dot product and proportion of composite spectrum is explained. 
#' Where 1 is perfect agreement and 0 is no agreement. If no spectral database
#' match has been made then the value is set to NA and this score will not be
#' used in calculating the average ranking.
#' \item \strong{"inSilico"} \emph{in silico} fragmentation data. Both the results of the 
#' \code{\link{metID.metFrag}} and \code{\link{metID.CFM}} functions. The total
#' proportion of the composite spectrum explained by each \emph{in silico}
#' fragmentation method (a value between 0-1) is used to rank the annotations.
#' If no \emph{in silico} fragmentation match has been made then a value of NA
#' is set and this score will not be used in calculating the average ranking.
#' \item \strong{"rtPred"} predicted retention time similarity. Annotations are ranked based on the 
#' retention time deviation from the predictive retention time model built using
#'  the function \code{\link{metID.rtPred}}. A relative score between 0-1 is
#'  calculated globally by taking the range of retention time deviation values.
#'  \item \strong{"chemSim"} chemical similarity score. The mean maximum 1st neighbour (connected by
#'  correlation \code{\link{metID.corrNetwork}} and/or spectral similarity
#'  \code{\link{metID.specSimNetwork}}) tanimoto chemical similarity scores
#'  calculated by \code{\link{metID.chemSim}} is used to rank annotations. 
#'  A relative score between 0-1 is calculated globally by taking the range of 
#'  mean maximum 1st neighbour chemical similarity scores.
#'  \item \strong{"pubMed"} crude literature based plausibility. The number of PMIDs returned by searching
#'  the compound name in PubMed. Number of returned PMIDs are used to generate
#'  a relative score ranking between 0-1. This aspect is highly reliant on the
#'  database name being the correct synonym to search the PubMed repository with.
#'  In an effort to ensure phospholipids are correctly search against PubMed
#'  a set of regular expressions has been created to identify common phospholipid
#'  annotations and use the compound class name rather than an abbreviation with 
#'  positional and fatty acid chain length information to obscure the number of
#'  pubmed abstract ids returned (see \code{\link{lipidAbbrev}}).
#'  
#'  This aspect is potentially time consuming (but only needs to be conducted once)
#'  as it complies closely with the NCBI recommendations from the section 
#'  \emph{"Frequency, Timing and Registration of E-utility URL Requests"} of book 
#'  \emph{"A General Introduction to the E-utilities"} by Eric Sayers \url{http://www.ncbi.nlm.nih.gov/books/NBK25497/}: 
#'  
#'  \strong{\emph{"In order not to overload the E-utility servers, NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure to comply with this policy may result in an IP address being blocked from accessing NCBI."}}
#'  
#'  This aspect is optional and will only work during these recommended times.
#'  However the function can optionally wait until the recommended time automatically.
#'  \item \strong{"substructure"} should the substructure score generated by the
#' \code{\link{dbProb}} function be used to rank possible annotations.
#' }
#' 
#' Depending on the availabilty of each of these pieces of evidence a mean 
#' annotation ranking score is calculated for every annotation and the best 
#' annotations can be automatically added. 
#' 
#' @source Sayers E. A General Introduction to the E-utilities. In: Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: http://www.ncbi.nlm.nih.gov/books/NBK25497
#' @export
setGeneric("metID.buildConsensus", function(object, ...) standardGeneric("metID.buildConsensus"))

setMethod("metID.buildConsensus", signature = "compMS2", function(object, 
include=NULL, metIDWeights=NULL, autoPossId=FALSE, minMeanBCscore=NULL, 
possContam=3, showPlot=TRUE, verbose=TRUE){
# error handling
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  if(length(DBanno(object)) == 0){
    stop('At minimum the functions metID.dbAnnotate and metID.dbProb must be run before building the annotation consensus')
  } 
  if(all(sapply(BestAnno(object), is.null))){
    stop('The function metID.dbProb must be run before building the annotation consensus.')  
  }
  # rbind all best annotations
  # index which bestAnnos have entries
  emptyIndx <- sapply(BestAnno(object), is.null) == FALSE
  if(all(emptyIndx == FALSE)){
    stop('none of composite spectra have any dbAnnotate matches')
  }
  # match bestAnno names to vertices
  # extract smiles codes 
  allBestAnnos <- do.call(rbind, BestAnno(object)[emptyIndx])
  nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
  allBestAnnos$specNamesTmp <- do.call(c, 
                                       mapply(rep, names(BestAnno(object)[emptyIndx]),
                                             each=nrowBestAnno))
  
  # remove any previous BC result columns
  allBestAnnos$BC_chemSim <- NULL
  allBestAnnos$BC_cfmScore <- NULL
  allBestAnnos$BC_massAccScore <- NULL
  allBestAnnos$BC_meanScore <- NULL
  allBestAnnos$BC_metFragScore <- NULL
  allBestAnnos$BC_nPubMedIds <- NULL
  allBestAnnos$BC_rtDevScore <- NULL
  allBestAnnos$BC_spectralDB <- NULL
  allBestAnnos$BC_substructure <- NULL
  allBestAnnos$BC_corrMolDesc <- NULL
  # inclusion names
  inclNames <- c('massAccuracy', 'spectralDB', 'inSilico', 'rtPred', 'chemSim', 
                 'pubMed', 'substructure', 'corrMolDesc')
  # scores to include
  if(is.null(include)){
    include <- inclNames
  }
  # check inclusion names are correct
  indxTmp <- include %in% inclNames 
  
  if(any(indxTmp == FALSE)){
  stop('The inclusion name(s):\n', paste0(include[indxTmp == FALSE], '\n'), 'are incorrect the inclusion names must be a combination of one or more of the following options:\n', 
       paste0(inclNames, '\n'))  
  }
  # check if metIDWeights is supplied that is matches length of include
  if(!is.null(metIDWeights)){
    if(length(metIDWeights) != length(include)){
      stop('metIDWeights argument must be the same length as the include argument')
    }
    if(sum(metIDWeights) >= 1.01 | sum(metIDWeights) <= 0.99){
      stop('The sum of the metIDWeights argument must be equal to 1') 
    }
  }
  # check to see if inclusion name appropriate based on CompMS2 content
  # 1. massAccuracy
  if(inclNames[1] %in% include){
  allBestAnnos$BC_massAccScore <- 1 - (abs(allBestAnnos$ppmMatch)/Parameters(object)$precursorPpm)
  }
  # 2. spectralDB
  if(inclNames[2] %in% include){
  if(length(object@spectralDB) == 0){
    stop('metID.matchSpectralDB must be run to calculate the spectralDB ranking score ("spectralDB")')
  } 
  # new score column
  allBestAnnos$BC_spectralDB <- 0
  dbMatchesTmp <- sapply(object@spectralDB, is.null) == F
  dbMatchIndx <- match(allBestAnnos$specNamesTmp, names(object@spectralDB)[dbMatchesTmp])
  # strsplit spectra
  for(i in unique(dbMatchIndx[!is.na(dbMatchIndx)])){
    dbSpectraDfTmp <- object@spectralDB[dbMatchesTmp][[i]]$dbSpectra
    dbSpecDfNames <- gsub('__.+', '', dbSpectraDfTmp[, 'compound_msp'])
    dbSpecDfNames <- strsplit(dbSpecDfNames, '')
    # namesTmp <- rep(1:length(dbSpecDfNames), sapply(dbSpecDfNames, length))
    # dbSpecDfNames <- unlist(dbSpecDfNames)
    # names(dbSpecDfNames) <- namesTmp
    indxTmp <- dbMatchIndx %in% i 
    allBestAnnos$BC_spectralDB[indxTmp] <- sapply(strsplit(allBestAnnos$DBname[indxTmp], ''), 
                                                  function(x){
      matchIndxTmp <- sapply(dbSpecDfNames, function(y) all(x %in% y))
      if(any(matchIndxTmp)){
      subSpecDfTmp <- dbSpectraDfTmp[matchIndxTmp, , drop=FALSE]
      maxDotProdIndx <- which.max(as.numeric(subSpecDfTmp[, 'dotProductScore']))
      dbMatchScore <- mean(as.numeric(subSpecDfTmp[maxDotProdIndx, c('dotProductScore', 'propTIC_explained')]))
      return(dbMatchScore)
      } else {
      return(0)  
      }})
  } # dbmatchindx for loop
  } # end spectralDB
  
  # 3. inSilico
  if(inclNames[3] %in% include){
  if(length(object@inSilico) == 0){
    stop('metID.metFrag and/or metID.CFM must be run to calculate the inSilico ranking score ("inSilico")')
  } 
  # i. metFrag 
  if('MetFrag' %in% names(object@inSilico)){
  # no result indx
  resultIndx <- which(sapply(object@inSilico$MetFrag, is.null) == FALSE)
  resultIndx <- resultIndx[sapply(object@inSilico$MetFrag[resultIndx], ncol) > 1]
  allMetFrag <- do.call(rbind, object@inSilico$MetFrag[resultIndx]) 
  nrowMetFrag <- sapply(object@inSilico$MetFrag[resultIndx], nrow)
  allMetFrag$compSpecName <- do.call(c, mapply(rep, 
                                               names(object@inSilico$MetFrag[resultIndx]), 
                                               each=nrowMetFrag))
  
  allMetFrag$propIntEx <- as.numeric(as.character(allMetFrag$propIntEx))
  allMetFrag$propIntEx[is.na(allMetFrag$propIntEx)] <- 0
  allMetFrag$Score <- as.numeric(as.character(allMetFrag$Score))
  allMetFrag$Score[is.na(allMetFrag$Score)] <- 0
  # average of prop intensity explained and metfrag score
  bcMetFragTmp <- apply(allMetFrag[, c('Score', 'propIntEx')], 1, mean)
  # add score to best annos
  allBestAnnos$BC_metFragScore <- 0
  indxTmp <- match(paste0(allBestAnnos$specNamesTmp, allBestAnnos$DBid),
                   paste0(allMetFrag$compSpecName, allMetFrag$DBid))
  # allBestAnnos$BC_metFragScore[!is.na(indxTmp)] <- allMetFrag$propIntEx[indxTmp[!is.na(indxTmp)]]
  # combined metfrag Score and prop of intensity explained
  allBestAnnos$BC_metFragScore[!is.na(indxTmp)] <- bcMetFragTmp[indxTmp[!is.na(indxTmp)]]
  }
  # ii. CFM 
  if('CFM' %in% names(object@inSilico)){
    # no result indx
    resultIndx <- which(sapply(object@inSilico$CFM, is.null) == FALSE)
    resultIndx <- resultIndx[sapply(object@inSilico$CFM[resultIndx], ncol) > 1]
    allCFM <- do.call(rbind, object@inSilico$CFM[resultIndx]) 
    nrowCFM <- sapply(object@inSilico$CFM[resultIndx], nrow)
    allCFM$compSpecName <- do.call(c, mapply(rep, 
                                                 names(object@inSilico$CFM[resultIndx]), 
                                                 each=nrowCFM))
    # remove duplicates
    allCFM <- allCFM[duplicated(paste0(allCFM$DBid, allCFM$compSpecName)) == FALSE, , drop=FALSE]
    # add score to best annos
    allBestAnnos$BC_cfmScore <- 0
    indxTmp <- match(paste0(allBestAnnos$specNamesTmp, allBestAnnos$DBid),
                     paste0(allCFM$compSpecName, allCFM$DBid))
    allBestAnnos$BC_cfmScore[!is.na(indxTmp)] <- allCFM$CFM_totPropEx[indxTmp[!is.na(indxTmp)]]
  }  
  } # end inSilico
  
  # 4. rtPred
  if(inclNames[4] %in% include){
    if(length(object@rtPred) == 0){
      stop('metID.rtPred must be run to calculate the predicted retention time ranking score ("rtPred")')
    }
    rtDevScore <- allBestAnnos$predRtDev
    indxTmp <- rtDevScore == ''
    maxRtDev <- max(abs(as.numeric(rtDevScore[indxTmp == FALSE])))
    rtDevScore[indxTmp] <- maxRtDev
    allBestAnnos$BC_rtDevScore <- 1 - (abs(as.numeric(rtDevScore))/maxRtDev)
  } # end rtPred

  # 5. chemSim
  if(inclNames[5] %in% include){
    if(('MMNNCSS' %in% colnames(allBestAnnos)) == FALSE){
      stop('metID.chemSim must be run to calculate the chemical similarity ranking score ("chemSim")')
    }
    allBestAnnos$MMNNCSS[is.na(allBestAnnos$MMNNCSS)] <- 0
    allBestAnnos$BC_chemSim <- allBestAnnos$MMNNCSS
    
    } # end chemSim
  
  # 6. PubMed
  if(inclNames[6] %in% include){
    if(is.null(allBestAnnos$nPubMedIds)){
    allBestAnnos$nPubMedIds <- 0
    uniqueDBnames <- unique(allBestAnnos$DBname)
    allLipidAbbrev <- uniqueDBnames
    
    for(i in 1:nrow(lipidAbbrev)){
      allLipidAbbrev[grepl(lipidAbbrev$regexpr[i], allLipidAbbrev, ignore.case = TRUE)] <- lipidAbbrev$Class[i]
    }
    # replace initial positional information which can mask the true number of PMIDs
    allLipidAbbrev <- gsub("^[0-9]+-|^[0-9,]+-|^[0-9,]+,-|^[0-9']+-|^[0-9]+-[A-z]-|^[A-z]-||^[0-9']+-[A-z]-|^[0-9]+[A-z]-|^[0-9]+:[0-9]+ |^[0-9]+:[0-9]+-|^[A-z][A-z]-", '', allLipidAbbrev)
    uniqueLipidAbbrev <- unique(allLipidAbbrev)
      
    recommendedTime <- TRUE
    eTz <- 'America/New_York'
    # eastern time and day of the week 
    eastTimeTmp <- as.POSIXct(format(Sys.time(), tz=eTz), tz=eTz)
    
    if(weekdays(eastTimeTmp) %in% c('Saturday', 'Sunday') == FALSE){
      startTime <- as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 21"), "%Y-%m-%d %H", tz=eTz)
      endTime <- as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 05"), "%Y-%m-%d %H", tz=eTz)
     # if not the weekend is it the right time to send entrez requests
     recommendedTime <- ifelse(eastTimeTmp <
                               as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 23:59:59"), "%Y-%m-%d %H:%M:%S", tz=eTz)
                               & eastTimeTmp > startTime
                               | eastTimeTmp < endTime
                               & eastTimeTmp >
                                 as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 00:00:01"), "%Y-%m-%d %H:%M:%S", tz=eTz),
                               T, F)
    }
    
    if(recommendedTime == FALSE){
      warning('\nIt is currently: ', weekdays(eastTimeTmp), ' ', eastTimeTmp, ' Eastern time\nThis is not the recommended time to send large requests to the NCBI servers:\nThe recommended time is on weekends or weekdays between the hours of 9 p.m. EST and 5 a.m. EST\n\nWould you like compMS2Miner to wait until the correct time and then proceed?\n\ntype (Y/N) and press [enter] to continue:', immediate. = TRUE)
    flush.console()
      waitTillCorrectTime <- readline()
      if(waitTillCorrectTime == 'Y'){
       while(recommendedTime == FALSE){
         
         timeRemaining <- difftime(startTime, eastTimeTmp, units = 'hours')
         # cat(paste("The time remaining for querying PubMed:", 
         #       formatC(as.numeric(timeRemaining),width=5,format='f',digits=2,flag='0'), 'hours. Please wait...\r'))
         # timeRemaining <- timeRemaining - 0.05
         # Sys.sleep(0.5)
         # timeRemaining <- as.numeric(round(timeRemaining, 2))
         # if sleep time is greater than time remaining then reduce the sleeptime
         sleepTimeTmp <- ifelse(timeRemaining > 600, 600, timeRemaining)
         # eastern time and day of the week 
         eastTimeTmp <- as.POSIXct(format(Sys.time(), tz=eTz), tz=eTz)
         message('Your time zone: ', Sys.time(), '.\nCurrent Eastern time: ', eastTimeTmp, '.\n', round(as.numeric(timeRemaining, units='hours'), 2), ' hours before start time (', startTime,' Eastern time).\nThis message will update every ', round(sleepTimeTmp/60, 0),' minutes.\n\n')
         flush.console()
         
         Sys.sleep(sleepTimeTmp) 
         
         if(weekdays(eastTimeTmp) %in% c('Saturday', 'Sunday') == FALSE){
           # if not the weekend is it the right time to send entrez requests 
           recommendedTime <- ifelse(eastTimeTmp < 
                                       as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 23:59:59"), "%Y-%m-%d %H:%M:%S", tz=eTz) 
                                     & eastTimeTmp > startTime 
                                     | eastTimeTmp < endTime 
                                     & eastTimeTmp >
                                       as.POSIXct(paste0(as.Date(eastTimeTmp, tz=eTz), " 00:00:01"), "%Y-%m-%d %H:%M:%S", tz=eTz),
                                     T, F)
         } 
        } # while loop recommended time
       } # if wait recommended time
      } # if recommended time is false
    
    # if recommended time is true
    if(recommendedTime == TRUE){
      # pmt <- proc.time() 
      # cut compound names in to 3 request then a second sleep
      # sysSleepTimes <- seq(3, length(uniqueLipidAbbrev), 3) 
      if(verbose == TRUE){ pb <- txtProgressBar(max=length(uniqueLipidAbbrev), style=3)}
      message('querying PubMed with ', prettyNum(length(uniqueLipidAbbrev), big.mark=','), ' database names, this process is time consuming but only needs to be carried out once, querying PubMed is limited to 3 queries per second...please wait.\n')
      flush.console()
      # test <- proc.time()
      nAbstractsPubMed <- vector('numeric', length(uniqueLipidAbbrev))
      for(i in 1:10){
        pmt1 <- proc.time()
        tmpDBname <- uniqueLipidAbbrev[i]
        if(verbose==TRUE){setTxtProgressBar(pb, i)}
         nAbsTmp <- PMIDsearch(tmpDBname, n=0)
         nAbstractsPubMed[i] <- ifelse(length(nAbsTmp) > 0, as.numeric(nAbsTmp), NA)
         pmt2 <- proc.time() - pmt1
         # sleep for whatever is left of the 1/3rd of a second
         slTime <- 0.33 - pmt2[3]
         Sys.sleep(ifelse(sign(slTime) == -1, 0, slTime))
         # cat(sign(slTime), '\n')
        # system sleep every 3 request for 2 seconds
        # if(i %in% sysSleepTimes){
        # Sys.sleep(1.5) 
        # }
      } # db names loop
       # proc.time() - test
    indxTmp <- match(allLipidAbbrev, uniqueLipidAbbrev) 
    nAbstractsPubMed <- nAbstractsPubMed[indxTmp]
    indxTmp <-  match(allBestAnnos$DBname, uniqueDBnames)
    allBestAnnos$nPubMedIds <- nAbstractsPubMed[indxTmp]
    } # if recommended time is true
    } # if nPMIDs column already
    nPubMedIdsTmp <- log(allBestAnnos$nPubMedIds)
    nPubMedIdsTmp[nPubMedIdsTmp == -Inf] <- 0
    allBestAnnos$BC_nPubMedIds <- nPubMedIdsTmp/max(nPubMedIdsTmp)
    if(all(is.na(allBestAnnos$BC_nPubMedIds))){
    allBestAnnos$nPubMedIds <- NULL
    allBestAnnos$BC_nPubMedIds <- NULL
    }
    } # end pubMed
  # 7. substructure score
  if(inclNames[7] %in% include){
    # allBestAnnos <- allBestAnnos[order(allBestAnnos$specNamesTmp), , drop=FALSE]
    # scaledScore <- tapply(allBestAnnos$subStrScore, allBestAnnos$specNamesTmp,
    #                       function(x){ifelse(x == 0, 0, x/max(x))})
    # scaledScore <- scaledScore[order(names(scaledScore))]
    # allBestAnnos$BC_substructure <- do.call(c, scaledScore)
    allBestAnnos$BC_substructure <- allBestAnnos$subStrScore
  } 
  
  # 8. correlation molecular descriptors internal databases
  if(inclNames[8] %in% include){
    allBestAnnos$BC_corrMolDesc <- apply(allBestAnnos[, grep('nCorrMolDesc', colnames(allBestAnnos))],
                                         1, sum)
    allBestAnnos <- allBestAnnos[order(allBestAnnos$specNamesTmp), , drop=FALSE]
    scaledScore <- tapply(allBestAnnos$BC_corrMolDesc, allBestAnnos$specNamesTmp,
                          function(x){ifelse(x == 0, 0, x/max(x))})
    scaledScore <- scaledScore[order(names(scaledScore))]
    allBestAnnos$BC_corrMolDesc <- do.call(c, scaledScore)
    # allBestAnnos$BC_corrMolDesc <- ifelse(allBestAnnos$BC_corrMolDesc > 0, 1, 0)
  }
  # mean annotation score
  if(is.null(metIDWeights)){
  allBestAnnos$BC_meanScore <- apply(allBestAnnos[, grep('^BC_', colnames(allBestAnnos)), drop=FALSE], 1, function(x){
    mean(as.numeric(x[!is.na(x) & x != '']))})
  } else { # weighted mean
    BCscores <- allBestAnnos[, grep('^BC_', colnames(allBestAnnos)), drop=FALSE]
    
    if(inclNames[1] %in% include){
    colnames(BCscores)[colnames(BCscores) == 'BC_massAccScore'] <- inclNames[1] 
    }
    if(inclNames[2] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_spectralDB'] <- inclNames[2] 
    }
    # if in silico then calc mean
    if(inclNames[3] %in% include){
    BCscores$inSilico  <- apply(BCscores[, grep('BC_metFragScore|BC_cfmScore', colnames(BCscores)), drop=FALSE], 1, mean)  
    BCscores$BC_metFragScore <- NULL
    BCscores$BC_cfmScore <- NULL
    }
    if(inclNames[4] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_rtDevScore'] <- inclNames[4] 
    }
    if(inclNames[5] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_chemSim'] <- inclNames[5] 
    }
    if(inclNames[6] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_nPubMedIds'] <- inclNames[6] 
    }
    if(inclNames[7] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_substructure'] <- inclNames[7] 
    }
    if(inclNames[8] %in% include){
      colnames(BCscores)[colnames(BCscores) == 'BC_corrMolDesc'] <- inclNames[8] 
    } 
    BCscores <- BCscores[, include, drop=FALSE]
    allBestAnnos$BC_meanScore <- apply(BCscores, 1, function(x) weighted.mean(x[!is.na(x) & x != ''], w=metIDWeights))   }
  
  allBestAnnos <- allBestAnnos[order(allBestAnnos$BC_meanScore, decreasing=TRUE), , drop=FALSE]
  # feedback plot BC_meanScore and interquartile ranges
  meanMeanBCscoreTmp <- mean(as.numeric(allBestAnnos$BC_meanScore))
  IQRsTmp <- c(meanMeanBCscoreTmp + IQR(allBestAnnos$BC_meanScore),
               meanMeanBCscoreTmp - IQR(allBestAnnos$BC_meanScore))
  
  if(showPlot == TRUE){
  plot(sort(allBestAnnos$BC_meanScore), xlab='Annotation rank', 
       ylab='mean consensus score', main='mean consensus scores')
  abline(h=IQRsTmp, col='red')
  abline(h=meanMeanBCscoreTmp, col='black')
  text(10, meanMeanBCscoreTmp, paste0('mean (', round(meanMeanBCscoreTmp, 2), ')'), pos=3)
  text(c(10, 10), IQRsTmp, c(paste0("+ IQR (", round(IQRsTmp[1], 2), ')'), paste0("- IQR (", round(IQRsTmp[2], 2), ')')), col="red", pos=3) 
  }
  message(sum(allBestAnnos$BC_meanScore >= ifelse(is.null(minMeanBCscore), IQRsTmp[1], minMeanBCscore)), ' annotations (n=', 
          length(unique(allBestAnnos$specNamesTmp[allBestAnnos$BC_meanScore >= ifelse(is.null(minMeanBCscore), IQRsTmp[1], minMeanBCscore)])),
          ' composite spectra) above ', 
          ifelse(is.null(minMeanBCscore), paste0('upper inter-quartile range (',
                 round(IQRsTmp[1], 2), ').\n'), 
                 paste0('minimum mean consensus score (', round(minMeanBCscore, 2),
                        ').\n')))
  flush.console()
  # automatic metID comment annotation
  if(autoPossId ==  TRUE){
    # message('CAUTION: do you wish to automatically add the most probable database annotations to the comments table (column possible_identity)?\n type (Y/N) and press [enter] to continue:')
    # flush.console()
    # overWritePossId <- readline()
    # if(overWritePossId == 'Y'){
      autoBestAnnosDf <- allBestAnnos[allBestAnnos$BC_meanScore >= ifelse(is.null(minMeanBCscore), IQRsTmp[1], minMeanBCscore), , drop=FALSE]
      # remove duplicates and keep top annotations
      autoBestAnnosDf <- autoBestAnnosDf[duplicated(autoBestAnnosDf$specNamesTmp) == FALSE, , drop=FALSE]
      # add possible annotations to metIDcomments
      alreadyAnno <- Comments(object)$compSpectrum[grepl('metID\\.matchSpectralDB', Comments(object)$user_comments)]
      if(length(alreadyAnno) > 0){
        message(length(alreadyAnno), ' were already identified by the metID.matchSpectralDB function.\n')
        flush.console()
        autoBestAnnosDf <- autoBestAnnosDf[{autoBestAnnosDf$specNamesTmp %in% alreadyAnno} == FALSE, , drop=FALSE]
      }
      if(nrow(autoBestAnnosDf) > 0){
        # met Id comments table
        metIDcomments <- Comments(object)
        
        # id possible contaminants
        possContaminants <- table(autoBestAnnosDf$DBname)
        if(any(possContaminants > possContam)){
          indxTmp <- possContaminants > possContam
          message('The following automatic possible annotations have been identified as possible contaminants (i.e. occuring more than ', possContam, ' times see possContam argument) and will be named "possible_contaminant" in the Comments table and flagged "metID.buildConsensus":\n', paste0(1:sum(indxTmp), '. ', names(possContaminants)[indxTmp], ' (occurs n=', possContaminants[indxTmp], ' times)\n'))
          flush.console()
          contamIndx <- (autoBestAnnosDf$DBname %in% names(possContaminants[indxTmp]))
          contamSpec <- autoBestAnnosDf$specNamesTmp[contamIndx]
          # add to comments
          indxTmp <- metIDcomments$compSpectrum %in% contamSpec
          metIDcomments$user_comments[indxTmp] <- 'metID.buildConsensus'
          metIDcomments$possible_identity[indxTmp] <- 'possible_contaminant'
          autoBestAnnosDf <- autoBestAnnosDf[contamIndx == FALSE, , drop=FALSE]
        }
        
        message(nrow(autoBestAnnosDf), ' composite spectra out of ', nrow(metIDcomments), ' total (', round((nrow(autoBestAnnosDf)/nrow(metIDcomments)) * 100, 2), '%) have been identified based on a mean consensus annotation score (>= ', round(ifelse(is.null(minMeanBCscore), IQRsTmp[1], minMeanBCscore), 2), ').\nThe following metrics were used to calculate the mean consensus annotation score:\n', paste0(include, collapse = '\n'), '\n')
        flush.console()
        message('These annotations will now be added to the "metID comments" table which can viewed in compMS2Explorer or accessed using the Comments function.\n')
        flush.console()
        # add possible annotations to metIDcomments
        indxTmp <- match(metIDcomments$compSpectrum, autoBestAnnosDf$specNamesTmp)
        metIDcomments$possible_identity[!is.na(indxTmp)] <- autoBestAnnosDf$DBname[indxTmp[!is.na(indxTmp)]] 
        metIDcomments$ESI_type[!is.na(indxTmp)] <- autoBestAnnosDf$ESI_type[indxTmp[!is.na(indxTmp)]]
        metIDcomments$user_comments[!is.na(indxTmp)] <- 'metID.buildConsensus'
        Comments(object) <- metIDcomments
      } else {
        warning('no mean consensus annotation scores were above the minimum score of ', round(ifelse(is.null(minMeanBCscore), IQRsTmp[1], minMeanBCscore), 2), ' consider reducing this parameter.\n')  
      }
  } # end autoPossId
  
  # split and add back to bestAnnos
  splitTmp <- allBestAnnos$specNamesTmp
  allBestAnnos$specNamesTmp <- NULL
  allBestAnnosList <- split(allBestAnnos, f=splitTmp)
  # add results back to object
  BestAnno(object) <- allBestAnnosList
  return(object)
}) # end function    
