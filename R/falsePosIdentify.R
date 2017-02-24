#' Identify and/or remove false positives from metID comments table based on a prohibited
#' keyword search
#' @param object a "compMS2" class object.
#' @param prohibKeyWords character a regular expression of prohibited keywords.
#' If any of these are discovered in any of the abstracts returned from the
#' PubMed Entrez search then they will be removed from the metID comments table.
#' @param n integer number of pubmedids to return per metID comments possible_identity.
#' (default = 50). Larger numbers may improve accuracy but increase computation time.
#' The maximum allowed value is 500 (limited by entrez system).
#' @param meanFreqPerAbs numeric minimum mean frequency of the summed frequencies
#' of the prohibited keywords. Default = 0.6 that is a mean summed frequency of 
#' any of the key words in the regular expression of 0.6 for it to be
#' considered. This limits the accidental removals of a true positive which
#' has a low word count frequency by chance.
#' @param removeFP logical if TRUE remove possible false positive annotations
#' from the possible_identity column of the metID comments table. If FALSE (default) the
#' possible false positive annotation names will remain but they will be
#' flagged in the user_comments column of the metID comments table as 
#' "possible false positive (falsePosRemoval)".
#' @param maxChar numeric maximum number of characters in cleaned abstract words to return.
#' @param verbose logical if TRUE display progress bars.
#' @export
setGeneric("falsePosIdentify", function(object, ...) standardGeneric("falsePosIdentify"))

setMethod("falsePosIdentify", signature = "compMS2", function(object, prohibKeyWords=NULL, 
                                                             n=50, 
                                                             meanFreqPerAbs=0.6,
                                                             removeFP=FALSE,
                                                             maxChar=50, 
                                                             verbose=TRUE, ...){
  if(!require(XML)){
    stop('The package XML is required to use this function')
  }
  if(!require(PubMedWordcloud)){
    stop('The package PubMedWordcloud is required to use this function')
  }
  metIDcomments <- Comments(object) 
  possIdIdx <- which(metIDcomments$possible_identity != '' & metIDcomments$possible_identity != 'no_annotations')
  uniqueDBnames <- unique(metIDcomments$possible_identity[possIdIdx])
  allLipidAbbrev <- uniqueDBnames
  
  for(i in 1:nrow(lipidAbbrev)){
    allLipidAbbrev[grepl(lipidAbbrev$regexpr[i], allLipidAbbrev, ignore.case = TRUE)] <- lipidAbbrev$Class[i]
  }
  # replace initial positional information which can mask the true number of PMIDs
  allLipidAbbrev <- gsub("^[0-9]+-|^[0-9,]+-|^[0-9,]+,-|^[0-9']+-|^[0-9]+-[A-z]-|^[A-z]-||^[0-9']+-[A-z]-|^[0-9]+[A-z]-|^[0-9]+:[0-9]+ |^[0-9]+:[0-9]+-|^[A-z][A-z]-| \\([A-z].+", '', allLipidAbbrev)
  uniqueLipidAbbrev <- unique(allLipidAbbrev)
  
  
  # alter lipids to their abbreviations
  prohibIdx <- vector('logical', length = length(uniqueLipidAbbrev))
  
  cat(paste0('Searching PubMed abstracts with the following prohibited keyword regular expression:\n"', prohibKeyWords, '"\n', length(uniqueLipidAbbrev), ' possible annotations from the metID comments table will be considered.\nAny annotation containing more than a mean summed frequency per abstract of ', 
             meanFreqPerAbs, ' of these prohibited keywords will be flagged in the user_comments column and/or optionally removed.\n'))
  if(verbose == TRUE){ pb <- txtProgressBar(max=length(uniqueLipidAbbrev), style=3)}
  # testing 
  maxFreq <- vector('numeric', length(uniqueLipidAbbrev))
  names(maxFreq) <- uniqueLipidAbbrev
  maxFreqNumAbs <- maxFreq
  for(i in 1:length(uniqueLipidAbbrev)){
  if(verbose==TRUE){setTxtProgressBar(pb, i)}
  pmidsTmp <- PMIDsearch(uniqueLipidAbbrev[i], n=n)
  if(length(pmidsTmp) > 1){
    absTmp <- getAbs(pmidsTmp[-1])
    if(length(absTmp) > 0){
      # clean abstracts
      ClAbs <- PubMedWordcloud::cleanAbstracts(absTmp)
      # only keep word which are less than max characters
      ClAbs <- ClAbs[which(sapply(as.character(ClAbs$word), nchar) < maxChar), , drop=FALSE]
      # remove compound name as usually the max word freq
      ClAbs <- ClAbs[grepl(uniqueLipidAbbrev[i], ClAbs$word, ignore.case = TRUE) == FALSE, , drop=FALSE]
      # ClAbs$freqPerc <- {ClAbs$freq/max(ClAbs$freq)} * 100
      # word present?
      wordPresIdx <- grepl(prohibKeyWords, ClAbs$word)
      if(any(wordPresIdx)){
        maxFreqNumAbs[i] <- sum(ClAbs$freq[wordPresIdx])/length(absTmp)
        prohibIdx[i] <- maxFreqNumAbs[i] >= meanFreqPerAbs
      }
     }
   }
  }
  cat('\ndone\n')
  
  if(any(prohibIdx)){
    cat(paste0(sum(prohibIdx), ' possible false positive annotations were detected using the following regular expression:\n', prohibKeyWords,'\n'))
    cat(paste0(1:sum(prohibIdx), '. ', names(maxFreqNumAbs[prohibIdx]), ' (mean sum freq/abs = ', round(maxFreqNumAbs[prohibIdx], 2), ')\n'))
    remAnnos <- uniqueDBnames[allLipidAbbrev %in% uniqueLipidAbbrev[prohibIdx]]
    remIdx <- metIDcomments$possible_identity %in% remAnnos
    # remove all previous false positive flags
    metIDcomments$user_comments <- gsub('possible false positive (falsePosRemoval) ', '', metIDcomments$user_comments)
    metIDcomments$user_comments[remIdx] <- paste0('possible false positive (falsePosRemoval) ', metIDcomments$user_comments[remIdx])
    if(removeFP == TRUE){
      metIDcomments$user_comments[remIdx] <- 'false positive removed by falsePosRemoval'
      metIDcomments$possible_identity[remIdx] <- ''
    }
  }
  Comments(object) <- metIDcomments
  return(object)
}) # end function
