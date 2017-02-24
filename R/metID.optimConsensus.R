#' differential evolution build consensus weight optimization
#' 
#' @details uses the package \code{\link{DEoptim}} to calculate the optimum
#' weighting of all included consensus scores to accurately rank the known
#' annotations (taken from the "metID comments" table in \code{\link{compMS2Explorer}}). These global parameters can then be used to rank the unknown annotations of other
#' unannotated composite spectra based on the optimum weighted mean consensus score.
#' The annotations above a certain score (minMeanBCscore) can also be automatically
#' added to the "metID comments" table. The ongoing differential evolution process will appear in a plot window with a loess fit line in red highlighting any reduction
#' in the mean rank of the training set annotations (from "metID comments" table) as
#' the genetic process evolves. This metaheuristic global optimization process
#' can help to maximise the parameters for accurate metabolite annotation and ranking.
#' 
#' @param object a "compMS2" class object.
#' @param include character vector of 7 options to build consensus combinatorial
#' metabolite identification see Details below for a description of each. If
#' specific options are not supplied as a character vector then the default 
#' is to consider all 7. i.e. 
#' c('massAccuracy', 'spectralDB', 'inSilico', 'rtPred', 'chemSim', 'pubMed', 'substructure').
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
#' @param specDbOnly logical if TRUE then only spectra matched to a spectral
#' database will be considered. These annotations are identified by the flag
#' "metID.matchSpectralDB" in the comments table. The default FALSE means that
#' all metabolites in the metID comments table will be considered.
#' @param  popSize numeric number of population members (see the NP argument in \code{\link{DEoptim.control}}). The default is 10 * length(include).
#' @param itermax numeric the maximum number of iterations (population generation)
#' allowed. (default = 100). See \code{\link{DEoptim.control}} for further details.
#' @param plotInterval the number of iterations before plotting the algorithms progress (default=20). Smaller values may slightly slow the process.
#' @param ... additional arguments to \code{\link{DEoptim.control}}.
#' @examples 
#' compMS2Example <- metID(compMS2Example, 'optimConsensus')
#' @references 
#' \enumerate{
#' \item David Ardia, Katharine M. Mullen, Brian G. Peterson, Joshua Ulrich (2015).
#' 'DEoptim': Differential Evolution in 'R'. version 2.2-3.
#' 
#' \item Katharine Mullen, David Ardia, David Gil, Donald Windover, James Cline
#' (2011). 'DEoptim': An R Package for Global Optimization by Differential
#' Evolution. Journal of Statistical Software, 40(6), 1-26. URL
#' http://www.jstatsoft.org/v40/i06/.
#' 
#' \item Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G. (2010).
#' Differential Evolution with 'DEoptim': An Application to Non-Convex Portfolio
#' Optimization. The R Journal, 3(1), 27-34. URL
#' http://journal.r-project.org/archive/2011-1/2011-1_index.html.
#' 
#' \item Ardia, D., Ospina Arango, N., Giraldo Gomez, N. (2010). Jump-Diffusion
#' Calibration using Differential Evolution. Wilmott Magazine, Issue 55
#' (September), 76-79. URL http://www.wilmott.com/.
#' 
#'  \item Kenneth V. Price, Rainer M. Storn and Jouni A. Lampinen (2006). 
#'  Differential Evolution - A Practical Approach to Global Optimization. Berlin Heidelberg:
#'   Springer-Verlag. ISBN 3540209506. 
#'   }
#' @export
setGeneric("metID.optimConsensus", function(object, ...) standardGeneric("metID.optimConsensus"))

setMethod("metID.optimConsensus", signature = "compMS2", function(object, include=NULL, 
autoPossId=FALSE, minMeanBCscore=NULL, possContam=3, specDbOnly=FALSE, popSize=NULL, itermax=100, 
plotInterval=20, ...){

  # error handling
  if(!require(DEoptim)){
    stop('package DEoptim is required to utilize this function')
  }
  # if(!require(RcppDE)){
  #   stop('package RcppDE is required to utilize this function')
  # }
  # inclusion names
  inclNames <- c('massAccuracy', 'spectralDB', 'inSilico', 'rtPred', 'chemSim',
                 'pubMed', 'substructure', 'corrMolDesc')
  # scores to include
  if(is.null(include)){
    include <- inclNames
  }
  message('Prior to score optimization.\n')
  flush.console()
  # run first buildConsensus
  object <- metID.buildConsensus(object, include=include, minMeanBCscore=minMeanBCscore,
                                 showPlot=FALSE)
  
  message('Starting global metabolite annotation consensus optimization.\n')
  flush.console()
  
  # index which bestAnnos have entries
  emptyIndx <- sapply(BestAnno(object), nrow) > 0
  if(all(emptyIndx == FALSE)){
    stop('none of the composite spectra have any dbAnnotate matches')
  }
  # best annos table
  bestAnnoDf <- do.call(rbind, BestAnno(object)[emptyIndx])
  nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
  bestAnnoDf$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[emptyIndx]), 
                                               each=nrowBestAnno))
  
  
  if(nrow(Comments(object)) == 0){
    stop('No comments have yet been made and therefore no possible_identity have been entered in the column.\n\nEither run the automatic chemical similarity metabolite identification (see ?metID.chemSim) or visualize with ?compMS2Explorer and add putative annotation names to the possible_identity column in the "met ID comments" tabbed panel') 
  }
  # extract comments table
  metIDcomments <- Comments(object)
  # match possible_identity to best anno table
  if(specDbOnly == TRUE){
    specDbIdx <- grepl('metID\\.matchSpectralDB', metIDcomments$user_comments)
    metIDcomments$possible_identity[specDbIdx == FALSE] <- '' 
  }
  indxTmp <- match(bestAnnoDf$specNamesTmp, metIDcomments$compSpectrum)
  bestAnnoDf$possible_identity <- metIDcomments$possible_identity[indxTmp]
  # which dbNames and possible_identity names match
  bestAnnoDf$nameMatch <- bestAnnoDf$DBname == bestAnnoDf$possible_identity 
  bestAnnoDf <- bestAnnoDf[order(bestAnnoDf$nameMatch, decreasing=TRUE), ]
  bestAnnoDf$nameMatch <- bestAnnoDf$nameMatch & duplicated(bestAnnoDf$specNamesTmp) == F & bestAnnoDf$DBname != ''
  if(all(bestAnnoDf$nameMatch == FALSE) & all(bestAnnoDf$possible_identity == '')){
    stop('No entries in the possible_identity column of the "Met ID comments" table.\n\nEither run the automatic chemical similarity metabolite identification (see ?metID.chemSim) or visualize with ?compMS2Explorer and add putative annotation names to the possible_identity column in the "met ID comments" tabbed panel')
  }
  # subset to include only those compspec with name matches
  bestAnnoDf <- bestAnnoDf[bestAnnoDf$specNamesTmp %in% bestAnnoDf$specNamesTmp[bestAnnoDf$nameMatch], , drop=FALSE]
  # replace column names
  if(inclNames[1] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_massAccScore'] <- inclNames[1] 
  }
  if(inclNames[2] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_spectralDB'] <- inclNames[2] 
  }
  # if in silico then calc mean
  if(inclNames[3] %in% include){
    bestAnnoDf$inSilico  <- apply(bestAnnoDf[, grep('BC_metFragScore|BC_cfmScore', colnames(bestAnnoDf)), drop=FALSE], 1, mean)  
    bestAnnoDf$BC_metFragScore <- NULL
    bestAnnoDf$BC_cfmScore <- NULL
  }
  if(inclNames[4] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_rtDevScore'] <- inclNames[4] 
  }
  if(inclNames[5] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_chemSim'] <- inclNames[5] 
  }
  if(inclNames[6] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_nPubMedIds'] <- inclNames[6] 
  }
  if(inclNames[7] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_substructure'] <- inclNames[7] 
  }
  if(inclNames[8] %in% include){
    colnames(bestAnnoDf)[colnames(bestAnnoDf) == 'BC_corrMolDesc'] <- inclNames[8] 
  } 
  # lower and upper bounds of each parameter to search
  lowerBound <- rep(0, length(include))
  upperBound <- rep(100, length(include))
  
  if(is.null(popSize)){
    popSize <- length(lowerBound) * 10
  }
  meanRank <<- numeric()

  # establish evaluation function
  evalFunc <- function(scoreWeights){
    weightsIt <- scoreWeights/ sum(scoreWeights)
    bestAnnoDf$BC_meanScore <- apply(bestAnnoDf[, include, drop=FALSE], 1, 
                                     function(x) weighted.mean(as.numeric(x[!is.na(x)|x!='']), w=weightsIt))  
    # sort by mean score
    bestAnnoDf <- bestAnnoDf[order(bestAnnoDf$BC_meanScore, decreasing = TRUE), , drop=FALSE]
    # find ranking of annotation for each spectrum
    idRankTmp <- tapply(bestAnnoDf$nameMatch, 
                        bestAnnoDf$specNamesTmp, 
                        function(x){which(x)})
    meanRank <<- c(meanRank, mean(idRankTmp))
    if(length(meanRank) %in% plotIntervalTmp){
    iterIndxTmp <- 1:length(meanRank) %in% iterationSeq
    if(length(meanRank) > plotInterval){
      dev.off()
    }
    #loess fit
    seqMeanRank <- 1:length(meanRank)
    loess_fit <- suppressWarnings(loess(meanRank ~ seqMeanRank))
    loess_line <- suppressWarnings(predict(loess_fit))
    plot(meanRank, type='l', ylab='mean rank all annotations',
    xlab='iteration number', ylim=c(0, max(meanRank)), main=paste0('current loess fit mean rank: ', round(loess_line[length(loess_line)], 2), ' (lowest: ', round(min(meanRank), 2),' highest: ', round(max(meanRank), 2), ')'), sub=paste0('generation number: ', ifelse(sum(iterIndxTmp) + 1 > nItersTmp, nItersTmp, sum(iterIndxTmp) + 1), ' of ', nItersTmp,  ' (nMembers = ', popSizeTmp, ')'))
    lines(1:length(meanRank), loess_line, col = "red", lwd=3)
    legend('bottomleft', paste0('training set = ', length(idRankTmp), ' annotations'))
    }
    setTxtProgressBar(pb, length(meanRank))
    return(sum(idRankTmp))
  }
  # # monitoring function
  popSizeTmp <<- popSize
  nItersTmp <<- itermax
  iterationSeq <<- seq(popSizeTmp, popSizeTmp * nItersTmp, popSizeTmp)
  plotIntervalTmp <<- seq(0, popSizeTmp * nItersTmp, plotInterval)
  pb <<- txtProgressBar(max=popSizeTmp * nItersTmp, style=3)
  # pmt <- proc.time()
  DEmodel <- DEoptim::DEoptim(evalFunc, lowerBound, upperBound,
                              DEoptim::DEoptim.control(NP=popSize, itermax=itermax,
                                                       trace=FALSE, ...))
  # proc.time()-pmt
  # optimum consensus weight identified
  ocWeights <- DEmodel$optim$bestmem/sum(DEmodel$optim$bestmem)
  message('The optimum metID weights identifed were as follows:\n')
  flush.console()
  # add to parameters
  for(i in 1:length(include)){
    Parameters(object)[, paste0('OC_', include[i])] <- round(ocWeights[i], 5)
    message(include[i], ': ', Parameters(object)[, paste0('OC_', include[i])], '\n')
    flush.console()
  }
  message('These will be saved in the compMS2 "Parameters" slot. Access using function Parameters().\n\n')
  flush.console()
  
  message('Building final consensus model using divergent evolution optimized weights...\n')
  flush.console()
  # run last buildConsensus with optimum weights
  object <- metID.buildConsensus(object, include=include, autoPossId=autoPossId, 
                                 minMeanBCscore=minMeanBCscore, 
                                 metIDWeights=ocWeights, possContam=possContam,
                                 showPlot=TRUE)
  
  return(object)
}) # end function
