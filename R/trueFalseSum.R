#' summary of true positive and false negative assignments for each score metric
#' 
#' @details The weighted average takes in to account the number of possible candidates
#' by weighting the mean correct rank by the number of possible candidates. e.g. In this
#' way a rank of 1 out of 10 possible structures will be ranked lower than a rank of 3 out
#' of 1,000 possible structures. The latter is a more important/impressive reflection of 
#' ranking ability of than the former. If a more extensive chemical database is
#' utilized such as PubMed compound for example then a far larger number of
#' candidates may be considered than a much smaller curated database such as HMDB. 
#' @param object a "compMS2" class object
#' @param n integer the number of top annotations to consider (default = 5).
#' i.e. if the correct structure is found in the top 5 structures then it is
#' considered a true positive.
#' @param minSimScore numeric chemical similarity score (values between 0 and 1).
#' This score is used to identify the compound contained in the metID comments
#' table (default = 1). The is performed rather than by name as there can
#' be multiple hits of the same structure in the top annotations. To check
#' if the metric is in roughly the right ballpark highly similar structures to 
#' the correct annotations can be considered by setting the argument to a lower 
#' Tanimoto chemical similarity score (e.g. 0.9).
#' @param specDbOnly logical if TRUE (default) then only the spectral database
#' annotations are used to calculate the true positive and false negative rates. If FALSE
#' then all possible_identities in the metID comments table are used.
#' @param verbose logical if TRUE display progress bars.
#' 
#' @return a list summarizing the true positive and false negative outputs, a 
#' summary plot of the results and also weighted average rank for each score metric.
#' 
#' @export
setGeneric("trueFalseSum", function(object, ...) standardGeneric("trueFalseSum"))

setMethod("trueFalseSum", signature = "compMS2", function(object, 
                                                          n=5, minSimScore=1,
                                                          specDbOnly=TRUE,
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
  # comments table
  metIDcomments <- Comments(object)
  if(all(metIDcomments$possible_identity == '')){
    stop('There are no compound names in the metID comments table column "possible_identity".')
  }
  if(specDbOnly == TRUE){
  corrAnnoIdx <- grepl('metID\\.matchSpectralDB', metIDcomments$user_comments)
  if(all(corrAnnoIdx == FALSE)){
    stop('There are no annotations which were matched by the metID.matchSpectralDB function.\n', 
         'The metID comments table column "user_comments" must contain the flag "metID.matchSpectralDB" make sure these are present and have not been edited/deleted.')
  }
  corrAnno <- metIDcomments[corrAnnoIdx, , drop=FALSE]
  } else {
  corrAnno <- metIDcomments[metIDcomments$possible_identity != '', , drop=FALSE]
  }
  # unique identifier 
  corrAnno$uniqueId <- paste0(corrAnno$compSpectrum, '_', corrAnno$possible_identity)
  # extract best anno information
  emptyIndx <- sapply(BestAnno(object), nrow) > 0
  if(all(emptyIndx == FALSE)){
    stop('none of composite spectra have any dbAnnotate matches')
  }
  # match to correct anno
  nameIdx <- names(BestAnno(object)) %in% corrAnno$compSpectrum
  # extract the top n
  smilesDf <- do.call(rbind, BestAnno(object)[emptyIndx & nameIdx])
  smilesDf$specNamesTmp <- gsub('\\..+', '', row.names(smilesDf))
  smilesDf$uniqueId <- paste0(smilesDf$specNamesTmp, '_', smilesDf$DBname)
  bcCols <- grepl('^BC_', colnames(smilesDf))
  # check build consensus run
  if(all(bcCols == FALSE)){
    stop('The metID.buildConsensus or metID.optimConsensus functions must be run before using this function.\n')
  }
  # chem fingerprints
  if(!'chemFP' %in% colnames(smilesDf)){
    stop('The function metID.chemSim must be run before using this function.')
  }
  idxTmp <- match(corrAnno$uniqueId, smilesDf$uniqueId)
  chemFpCorrAnno <- smilesDf$chemFP[idxTmp]
  names(chemFpCorrAnno) <- smilesDf$uniqueId[idxTmp]
  # res list
  resList <- vector('list', sum(bcCols))
  names(resList) <- colnames(smilesDf)[bcCols]
  # rank list
  rankList <- vector('list', sum(bcCols))
  names(rankList) <- colnames(smilesDf)[bcCols]
  # n unique chemFP
  corrAnno$uniqueChemFP <- 0
  if(verbose == TRUE){ pb <- txtProgressBar(max=nrow(corrAnno), style=3)}
  for(i in 1:nrow(corrAnno)){
    if(verbose==TRUE){setTxtProgressBar(pb, i)}
    smilesSub <- smilesDf[smilesDf$specNamesTmp %in% corrAnno$compSpectrum[i], , drop=FALSE]
    # chemFp matrix
    bAChemFP <-  do.call(rbind, lapply(smilesSub$chemFP, function(x){
      indxTmp <- as.numeric(strsplit(x, ';')[[1]])
      return(ifelse(1:1024 %in% indxTmp, 1, 0))}))
    # calculate chemical similarity
    quMol <- vector('numeric', 1024)
    quMol[as.numeric(strsplit(chemFpCorrAnno[i], ';')[[1]])] <- 1
    smilesSub$simMol <- ChemmineR::fpSim(quMol, bAChemFP, sorted = FALSE)
    # take in to account unique ESI adduct types
    esiTypesNum <- as.numeric(factor(smilesSub$ESI_type))
    # number of unique chemical fingerprints
    smilesSub$simMolEsiType <- smilesSub$simMol + esiTypesNum 
    corrAnno$uniqueChemFP[i] <- length(unique(smilesSub$simMolEsiType))
    # order by score metric and take tp/fn rate
    for(j in 1:length(resList)){
    smilesSub <- smilesSub[order(smilesSub[, names(resList)[j]], decreasing = TRUE), , drop=FALSE]
    # remove duplicate similarity scores
    smilesTop <- smilesSub[duplicated(smilesSub$simMolEsiType) == FALSE, , drop=FALSE]
    smilesTop <- smilesTop[smilesTop[, names(resList)[j]] != 0, , drop=FALSE]
    
    if(nrow(smilesTop) > 0){
    corrMetabIdx <- smilesTop$simMol >= minSimScore & smilesTop$ESI_type %in% corrAnno$ESI_type[i]
    resTmp <- ifelse(any(corrMetabIdx), 'FN', 'zeroScore')
    rankTmp <- which(corrMetabIdx)
    smilesTop <- smilesTop[1:ifelse(nrow(smilesTop) < n, nrow(smilesTop), n), , drop=FALSE] 
    corrMetabIdx <- smilesTop$simMol >= minSimScore & smilesTop$ESI_type %in% corrAnno$ESI_type[i]
    resTmp <- ifelse(any(corrMetabIdx), 'TP', resTmp)
    } else {
    resTmp <- 'zeroScore'
    rankTmp <- integer()
    }
    # if no score then return half the number of unique fingerprints/esi adducts
    rankTmp <- ifelse(length(rankTmp) == 0, round(corrAnno$uniqueChemFP[i]*0.5, 0), 
                      min(rankTmp))
    # rankTmp <- ifelse(length(rankTmp) == 0, NA, min(rankTmp))
    names(resTmp) <- names(chemFpCorrAnno[i])
    names(rankTmp) <- names(chemFpCorrAnno[i])
    # calc chem sim
    resList[[j]] <- c(resList[[j]], resTmp)
    rankList[[j]] <- c(rankList[[j]], rankTmp)
    }
  }
  # grouped barplot
  if(!require(ggplot2)){
    warning('you must install ggplot2 to view the barplot', immediate. = TRUE)
  } else {
    resDf <- data.frame(stringsAsFactors = TRUE)
    for(k in 1:length(resList)){
      summTmp <- table(resList[[k]])
      if(!'zeroScore' %in% names(summTmp)){
        summTmp <- c(summTmp, zeroScore=0) 
      }
      if(!'FN' %in% names(summTmp)){
        summTmp <- c(FN=0, summTmp) 
      }
      if(!'TP' %in% names(summTmp)){
        summTmp <- c(TP=0, summTmp) 
      }
      resTmp <-  data.frame(scoreMetric=rep(names(resList)[k], 3), 
                            result=names(summTmp), n=as.numeric(summTmp), 
                            stringsAsFactors = TRUE)
      resDf <- rbind(resDf, resTmp)
    }
    # plot it 
    p <- ggplot(resDf, aes(x=scoreMetric, y=n, fill=result)) + 
      geom_bar(stat="identity", position="dodge") + 
      ggtitle(paste0('score metric TP/FP summary (top ', n, ' min chemical Sim=', round(minSimScore, 2), ')'  )) +
      geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)+
      ylab("count") +
      coord_cartesian(ylim = c(0, nrow(corrAnno))) +
      theme(text=element_text(size=20), axis.text.x = element_text(angle=45, hjust=1), axis.ticks = element_line(colour = "black"), 
            legend.key = element_rect(colour = "grey80"), panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey50"), 
            panel.grid.major = element_line(colour = "grey90", size = 0.2), 
            panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
            strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
    print(p)
    }
  # true positive and false negative rates and weighted means
  tpFnRates <- vector('list', length(resList))
  names(tpFnRates) <- names(resList)
  
  for(i in 1:length(resList)){
  # weighted means 
  rankTmp <- rankList[[i]]
  resTmp <- resList[[i]]
  zeroScoreIdx <- resTmp == 'zeroScore'
  wMeanTmp <- round(weighted.mean(rankTmp, corrAnno$uniqueChemFP, na.rm=TRUE), 2)
  meanTmp <- round(mean(rankTmp, na.rm=TRUE), 2)
  medianTmp <- round(median(rankTmp, na.rm=TRUE), 2)
  # TP FP rates
  resTmp <- resList[[i]] 
  allzIdx <- resTmp != 'zeroScore'
  resTmp <- resTmp[allzIdx]
  tpR <- round({sum(resTmp == 'TP')/length(resTmp)} * 100, 2)
  fnR <- round({sum(resTmp == 'FN')/length(resTmp)} * 100, 2)
  tpFnRates[[i]][[paste0('TPrate', ifelse(any(allzIdx == FALSE), ' (if non zero)', ''))]] <- tpR
  tpFnRates[[i]][[paste0('FNrate', ifelse(any(allzIdx == FALSE), ' (if non zero)', ''))]] <- fnR
  tpFnRates[[i]][['weightedMeanRank']] <- wMeanTmp
  tpFnRates[[i]][['meanRank']] <- meanTmp
  tpFnRates[[i]][['medianRank']] <- medianTmp
  tpFnRates[[i]]['nZeroScore'] <- sum(zeroScoreIdx)
  tpFnRates[[i]]['percZeroScore'] <- round({sum(zeroScoreIdx)/nrow(corrAnno)} * 100, 2)
  } 
  print(tpFnRates)
  return(list(tpTnResults=resList, ranks=rankList, tpFnRates=tpFnRates, nUniqueChemFP=corrAnno$uniqueChemFP))
}) # end function
