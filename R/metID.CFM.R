#' Wrapper function for Competitive Fragmentation modelling (CFM) \emph{in silico} fragmentation software
#' 
#' @description The function will automatically run as a parallel computation is the compMS2 object was created in parallel.
#' 
#' @param object a "compMS2" class object.
#' @param featureSubSet character vector of composite spectra names (e.g. CC_1, CC_2 etc.) otherwise the default is to perform CFM fragmentation on all composite spectra.
#' @param keepTempFiles logical default = FALSE fraggraph-gen .csv output file will
#' be created as temporary files otherwise if TRUE file will be retained in subdirectories named by composite spectrum name.
#' @param minPropTicEx numeric minimum mean total ion current explained (default = 0.9)
#' the candidate with the highest proportion of the total ion current explained
#'  above this minimum will be automatically added to the Comments
#' table. The argument autoPossId must also be set to TRUE.
#' @param autoPossId logical if TRUE the function will automatically add the name
#'  of the top annotation based on mean total ion current explained and metFrag score
#'   (default = FALSE). Caution if TRUE
#'  this will overwrite any existing possible_identities in the "metID comments"
#'  table. This functionality is intended as an automatic annotation identification tool prior to thorough examination of the data in \code{\link{compMS2Explorer}}.
#'  The intention is that automatic annotations can be used in the metID.rtPred
#'  retention prediction function as part of a seamless first-pass workflow.
#' @param possContam numeric how many times does a possible annotation have
#'  to appear in the automatically generated possible annotations for it to be
#'  considered a contaminant and therefore not added to the "metID comment" table (default = 3, i.e. if a database name appears more than 3 times in the 
#'  automatic annotation table it will be removed).
#' @param verbose logical if TRUE display progress bars. 
#' @source 
#' \enumerate{
#' \item Allen F, Pon A, Wilson M, Greiner R, and Wishart D. CFM-ID: a web server for annotation, spectrum prediction and metabolite identification from tandem mass spectra. Nucleic Acids Res. June 2014. \url{http://nar.oxfordjournals.org/content/early/2014/06/03/nar.gku436.full}.
#' \item fraggraph-gen.exe file in extdata downloaded (2016/07/09, cfm-id-2.2_win32.zip) from \url{https://sourceforge.net/p/cfm-id/wiki/Home/}. 
#' \item lpsolve.dll in extdata downloaded (2016/07/09, lp_solve_5.5.2.3_IDE_Setup.exe) from \url{https://sourceforge.net/projects/lpsolve/}.
#' }
#' @export
setGeneric("metID.CFM", function(object, ...) standardGeneric("metID.CFM"))

setMethod("metID.CFM", signature = "compMS2", function(object, 
                                                      featureSubSet = NULL,  
                                                      keepTempFiles = FALSE,
                                                      minPropTicEx=0.9,
                                                      autoPossId=FALSE,
                                                      possContam=3,
                                                      verbose=TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("BestAnnotations function has not yet been run or no best annotations have been selected")
  } else if (all(sapply(BestAnno(object), is.null))){
    stop("BestAnnotations function has not yet been run or no best annotations have been selected")
  } 
    # if (is.null(fragGraphGenExe)){
    #   tcltk::tkmessageBox(message = 'See "http://sourceforge.net/p/cfm-id/wiki/Home/" Competitive Fragmentation Modelling download zipfile and then finally locate fraggraph-gen.exe file (the file lpsolve55.dll must also be in the same directory as the fraggraph-gen.exe file)')
    #   
    #   fragGraphGenExe <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{Executable} {.exe}} {{All files} *}",
    #                                                    title="select your fraggraph-gen.exe file"))
    #   # paramsLog <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{log file} {.log}} {{All files} *}",
    #   #                                                   title="select the params_log file"))
    #   
    # }
  fragGraphGenExe <- system.file('extdata', 'fraggraph-gen.exe', 
                                 package='compMS2Miner')
    
    # if CFM file empty then create list
    if(length(CFM(object)) == 0){
      CFMResTmp <- vector("list", length(compSpectra(object)))
      names(CFMResTmp) <- names(compSpectra(object))
      CFM(object) <- CFMResTmp
    }
    
    if(is.null(featureSubSet)){
      featureSubSet <- names(compSpectra(object))[sapply(BestAnno(object), nrow) > 0]
    }
    featureSubSet.indx <- featureSubSet %in% names(compSpectra(object)) 
    if(any(featureSubSet.indx == FALSE)){
      stop("The following composite spectra subset names do not match :", 
           paste0(sapply(featsureSubSet[featureSubSet.indx == FALSE], message), "\n"))
    }
    
    
    #     if(verbose == TRUE){ pb <- txtProgressBar(min=0,max=length(featureSubSet),style=3)}
    featureSubSet <- which(names(compSpectra(object)) %in% featureSubSet)
    bestAnnoAll <- do.call(rbind, BestAnno(object)[featureSubSet])
    nrowBestAnno <- sapply(BestAnno(object)[featureSubSet], nrow)
    bestAnnoAll$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[featureSubSet]),
                                                  each=nrowBestAnno, SIMPLIFY=FALSE))
    bestAnnoSub <- bestAnnoAll[bestAnnoAll$SubStr_type == '', ] 
    bestAnnoSub <- bestAnnoSub[bestAnnoSub$ESI_type %in% ifelse(Parameters(object)$mode == 'pos', '[M+H]+', '[M-H]-'), ]
    # which query spectra still have smiles
    querySpecIndxTmp <- match(bestAnnoSub$specNamesTmp, names(compSpectra(object)))  
    bestAnnoSub <- bestAnnoSub[!is.na(querySpecIndxTmp), ]
    querySpecIndxTmp <- querySpecIndxTmp[!is.na(querySpecIndxTmp)]
    featureSubSet <- unique(querySpecIndxTmp)
    compSpecAll <- do.call(rbind, compSpectra(object)[featureSubSet])
    nrowBestAnno <- sapply(compSpectra(object)[featureSubSet], nrow)
    specNamesTmp <- do.call(c, mapply(rep, names(compSpectra(object)[featureSubSet]),
                                                  each=nrowBestAnno, SIMPLIFY=FALSE))
    compSpecAll <- compSpecAll[, c('mass', 'intensity')]
    compSpecAll <- cbind(compSpecAll, 
                         featureSubSet=rep(featureSubSet, table(specNamesTmp)))
    querySpecDbId <- tapply(querySpecIndxTmp, bestAnnoSub$DBid, paste0, collapse=' ')
    bestAnnoSub <- bestAnnoSub[duplicated(bestAnnoSub$DBid) == FALSE, ]
    bestAnnoSub <- cbind(bestAnnoSub, querySpecDbId=querySpecDbId[match(bestAnnoSub$DBid, names(querySpecDbId))])
    # create temporary directory
    
    message('calculating fragmentation graphs for ', nrow(bestAnnoSub), ' unique SMILES codes...please wait\n')
    flush.console()
    pmt <- proc.time()
    # parallel 
    if(Parameters(object)$nCores > 0){
      if(!require(foreach)){
        stop('package foreach must be installed to use this function in parallel')
      }
      if(!require(doSNOW)){
        stop('package doSNOW must be installed to use this function in parallel')
      }
      message(paste0("Starting SNOW cluster with ", Parameters(object)$nCores, " local sockets..."))
      flush.console()
      cl <- parallel::makeCluster(Parameters(object)$nCores, outfile='') 
      doSNOW::registerDoSNOW(cl)
      
      progress <- function(n){if(n %% 100 == 0){cat(paste0(n, ' of ', nrow(bestAnnoSub),
                                                          ' fragmentation graphs complete.\n'))}}
      
      if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
      
      # foreach and dopar from foreach package
      cfmFragResults <- foreach(i=1:nrow(bestAnnoSub), .options.snow=opts) %dopar% {
                cfmFragGraphGen(bestAnnoSubRow=bestAnnoSub[i, , drop=FALSE],
                                compSpecAll=compSpecAll,
                                keepTempFiles=keepTempFiles, 
                                fragGraphGenExe=fragGraphGenExe,
                                mode=Parameters(object)$mode,
                                frag_mzabs=Parameters(object)$Frag_mzabs)
      }
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      
    } else {  
    if(verbose == TRUE){ pb <- txtProgressBar(max=nrow(bestAnnoSub), style=3)}
    cfmFragResults <- vector('list', nrow(bestAnnoSub))
    for(i in 1:nrow(bestAnnoSub)){
      # progress bar
      if(verbose==TRUE){setTxtProgressBar(pb, i)}
      
      cfmFragResults[[i]] <- cfmFragGraphGen(bestAnnoSubRow=bestAnnoSub[i, , drop=FALSE],
                                             compSpecAll=compSpecAll,
                                             keepTempFiles=keepTempFiles, 
                                             fragGraphGenExe=fragGraphGenExe,
                                             mode=Parameters(object)$mode,
                                             frag_mzabs=Parameters(object)$Frag_mzabs)
    } # end loop
    } # end if parallel cond
    proc.time() - pmt
    cfmFragResults <- unlist(cfmFragResults, recursive = FALSE)
    cfmFragResults <- tapply(cfmFragResults, names(cfmFragResults), function(x){
      do.call(rbind, x)})
      # prev results
      prevResIndx <- sapply(object@inSilico$CFM[as.numeric(names(cfmFragResults))], is.data.frame)
      match(names(object@inSilico$CFM), names(cfmFragResults))
      if(any(prevResIndx == FALSE)){
        object@inSilico$CFM[as.numeric(names(cfmFragResults))[prevResIndx == FALSE]] <- cfmFragResults[prevResIndx == FALSE] 
      } 
      if(any(prevResIndx)){
        object@inSilico$CFM[as.numeric(names(cfmFragResults))[prevResIndx]] <- lapply(as.numeric(names(cfmFragResults))[prevResIndx], function(x) rbind(object@inSilico$CFM[[x]], cfmFragResults[[as.character(x)]])) 
      }
    proc.time() - pmt
    
    # if necessary then add to comments table
    if(autoPossId ==  TRUE){
      # message('CAUTION: do you wish to automatically add the most probable database annotations to the comments table (column possible_identity)?\n type (Y/N) and press [enter] to continue (NB. matchSpectralDB matches will not be overwritten):')
      # flush.console()
      # overWritePossId <- readline()
      # if(overWritePossId == 'Y'){
       bestCandidates <-  sapply(CFM(object), function(x){
          bestCand <- matrix('', ncol=2)
          if(!is.null(x)){
            if(ncol(x) > 1){
              propExTmp <- as.numeric(as.character(x[, 'CFM_totPropEx']))
              propExTmp[is.na(propExTmp)] <- 0
              
              if(any(propExTmp >= minPropTicEx)){
                bestCand <- x[which.max(propExTmp), 'DBname']
                bestCand <- cbind(bestCand, ifelse(Parameters(object)$mode == 'neg', 
                                                   '[M-H]-', '[M+H]+'))
              }
            }
          }
          return(bestCand)})
        
        if(any(bestCandidates[1, ] != '')){
          # met Id comments table
          metIDcomments <- Comments(object)
          bestCandidates <- bestCandidates[, bestCandidates[1, ] != '', drop=FALSE]
          # id possible contaminants
          possContaminants <- table(bestCandidates[1, ])
          if(any(possContaminants > possContam)){
            indxTmp <- possContaminants > possContam
            message('The following automatic possible annotations have been identified as possible contaminants (i.e. occuring more than ', possContam, ' times see possContam argument) and will be named "possible_contaminant" in the Comments table and flagged "metID.CFM":\n', paste0(1:sum(indxTmp), '. ', names(possContaminants)[indxTmp], ' (occurs n=', possContaminants[indxTmp], ' times)\n'))
            flush.console()
            contamIndx <- (bestCandidates[1, ] %in% names(possContaminants[indxTmp]))
            contamSpec <- colnames(bestCandidates)[contamIndx]
            # add to comments
            indxTmp <- metIDcomments$compSpectrum %in% contamSpec
            metIDcomments$user_comments[indxTmp] <- 'metID.CFM'
            metIDcomments$possible_identity[indxTmp] <- 'possible_contaminant'
            bestCandidates <- bestCandidates[, contamIndx == FALSE, drop=FALSE]
          }
          
          message(ncol(bestCandidates), ' composite spectra out of ', nrow(metIDcomments), ' total (', round((ncol(bestCandidates)/nrow(metIDcomments)) * 100, 2), '%) have been annotated based on a minimum proportion of the total ion current explained (>= ', round(minPropTicEx, 2), ').\n\n')
          flush.console()
          message('These annotations will now be added to the "metID comments" table in compMS2Explorer and flagged as "metID.CFM".\n')
          flush.console()
          # add possible annotations to metIDcomments
          alreadyAnno <- Comments(object)$compSpectrum[grepl('metID\\.matchSpectralDB', Comments(object)$user_comments)]
          if(length(alreadyAnno) > 0){
            bestCandidates <- bestCandidates[, {colnames(bestCandidates) %in% alreadyAnno} == FALSE, drop=FALSE]
          }
          if(ncol(bestCandidates) > 0){
          indxTmp <- match(metIDcomments$compSpectrum, colnames(bestCandidates))
          metIDcomments$possible_identity[!is.na(indxTmp)] <- bestCandidates[1, indxTmp[!is.na(indxTmp)]]
          metIDcomments$ESI_type[!is.na(indxTmp)] <- bestCandidates[2, indxTmp[!is.na(indxTmp)]]
          metIDcomments$user_comments[!is.na(indxTmp)] <- 'metID.CFM'
          Comments(object) <- metIDcomments
        } else {
          message('no candidates were automatically annotated by metID.CFM try lowering the minimum proportion of the TIC explained parameter (minPropTicEx).\n')
          flush.console()
        }
        }
      }
    
    return(object)
}) # end function
