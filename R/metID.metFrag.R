#' metFrag \emph{in silico} fragmentation query wrapper function 
#' 
#' @description performs metFrag (msbi.ipb-halle.de/MetFrag/) \emph{in silico} combinatorial
#' fragmentation. Local chemical structure data files (.sdf) are created from 
#' most probable annotation canonical SMILES codes. Temporary local sdf files 
#' and metfrag parameter files (.txt) are  created on a composite spectrum by 
#' composite spectrum basis and
#' \emph{in silico} fragmentation performed. Results are read back into R and stored in the compMS2 class object as results tables. In addition the temporary sdf,
#' .txt and .csv metfrag results files may also optionally be kept (keepTempFiles = TRUE) and are saved in a subdirectory structure (see \emph{details}). The function will automatically run as a parallel computation is the compMS2 object was created in parallel. The adduct types and codes which MetFrag command line tool 
#' will work with are included in a default internal table. See data(metFragAdducts)
#'  the following electrospray adducts are contained within:
#' positive mode: '[M+H]+', '[M+NH4]+', '[M+Na]+', '[M+K]+'
#' negative mode: '[M-H]-', '[M+Cl]-', '[M-H+HCOOH]-', '[M-H+CH3COOH]-' 
#' All other database annotations of other electrospray adducts will be discarded.
#' If the format of the metFragAdducts table is correctly followed then additional
#' MetFrag adduct types can be added see View(metFragAdducts). This ensures that
#' the adduct types are customizable and can be modified to incorporate future
#' availability. See argument metFragAdducts below.
#'@param object a "compMS2" class object.
#' @param featureSubSet character vector of composite spectra names (e.g. CC_1, CC_2 etc.) otherwise the default is to perform metFrag queries on all composite spectra.
#' @param adductTable data.frame containing adduct names, MetFrag number codes
#' and ionization mode. A default data.frame is internal to the package but this
#' table is fully customizable if future adduct types are added. See View(metFragAdducts)
#' to see the format of this table. N.B. The metFrag adduct names must perfectly match
#' those supplied to the \code{\link{adduct2mass}} function during database
#' annotation using \code{\link{metID.dbAnnotate}}.
#' @param keepTempFiles logical default = FALSE, txt, sdf and csv results files will
#' be created as temporary files otherwise if TRUE files will be retained in named subdirectories (see details).
#' @param maxTreeDepth numeric maximum tree depth (default = 2). If 2 then
#' fragments of fragments are also considered but will increase computation time.
#' @param frag_mzabs numeric delta predicted-observed fragment mass accuracy for matching.
#' @param minMetFragScore numeric minimum mean total ion current explained and metFrag score (default = 0.9)
#' the candidate with the highest score above this minimum will be automatically added to the Comments
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
#'@details if keepTempFiles = FALSE, Results directories are generated in the current working directory: for each MS1 feature matched to MS2 data a results 
#'                      directory is created named ("MetFrag_results"). Subdirectories are then created within the results directory named after each composite spectrum name.
#'                      Assuming MetFrag returned any results the following files should appear within the subdirectories:
#' \enumerate{
#' \item MetFrag parameters files : MetFrag parameter files (.txt) are saved in each result directory. See \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/} for details.
#' \item localSdf file : the local sdf (chemical structure data files) are saved in 
#'                     each result directory. This file is used as a local database for MetFrag \emph{in silico} fragmentation.
#' \item result .csv file: the metFrag results are returned as a comma seperated values text file (.csv).
#' }
#' @return a compMS2 class object containing metFrag \emph{in silico} fragmentation 
#' results which can be visualized in \code{\link{compMS2Explorer}}.
#' @source \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/} developed based on the command line .jar file (MetFrag2.3-CL.jar) downloaded on 2016/07/12
#' \enumerate{
#' \item MetFrag relaunched: incorporating strategies beyond in silico fragmentation: C Ruttkies, E L Schymanski, S Wolf, J Hollender, S Neumann Journal of Cheminformatics 2016 8:3
#' \item In silico fragmentation for computer assisted identification of metabolite mass spectra: S Wolf, S Schmidt, M Müller-Hannemann, S Neumann BMC bioinformatics 11 (1), 148
#' }
#' @export
setGeneric("metID.metFrag", function(object, ...) standardGeneric("metID.metFrag"))

setMethod("metID.metFrag", signature = "compMS2", function(object,  
                                                          featureSubSet=NULL,
                                                          adductTable=metFragAdducts,
                                                          keepTempFiles=FALSE,
                                                          maxTreeDepth=2,
                                                          frag_mzabs=0.05, 
                                                          minMetFragScore=0.9,
                                                          autoPossId=FALSE,
                                                          possContam=3,
                                                          verbose=TRUE){ 
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("metID.dbProb function has not yet been run or no probable annotations have been selected")
  } else if (all(sapply(BestAnno(object), is.null))){
    stop("metID.dbProb function has not yet been run or no probable annotations have been selected")
  } else if (!require(ChemmineR)){
    stop("the ChemmineR package must be installed from the Bioconductor repository to proceed.")
  } else if (!require(ChemmineOB)){
    stop("the ChemmineOB package must be installed from the Bioconductor repository to proceed.")
  } 
    
  metFragJar <- system.file('extdata', 'MetFrag2.3-CL.jar', package='compMS2Miner')
  Parameters(object)$metFrag_frag_mzabs <- frag_mzabs
  Parameters(object)$metFrag_maxTreeDepth <- maxTreeDepth
    # if metFrag file empty then create list
    if(length(MetFrag(object)) == 0){
     metFragResTmp <- vector("list", length(compSpectra(object)))
     names(metFragResTmp) <- names(compSpectra(object))
     MetFrag(object) <- metFragResTmp
    }
    
    if(is.null(featureSubSet)){
      featureSubSet <- names(compSpectra(object))[sapply(BestAnno(object), nrow) > 0]
    }
    featureSubSet.indx <- featureSubSet %in% names(compSpectra(object)) 
    if(any(featureSubSet.indx == FALSE)){
      stop("The following composite spectra subset names do not match :", 
           paste0(sapply(featsureSubSet[featureSubSet.indx == FALSE], message), "\n"))
    }
    
    bestAnnoAll <- do.call(rbind, BestAnno(object)[featureSubSet])
    nrowBestAnno <- sapply(BestAnno(object)[featureSubSet], nrow)
    bestAnnoAll$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[featureSubSet]),
                                               each=nrowBestAnno, SIMPLIFY=FALSE))
    
    bestAnnoSub <- bestAnnoAll[bestAnnoAll$SubStr_type == '', ]
    # if ion list for each mode
    if(Parameters(object)$mode == 'pos'){
    esiList <- adductTable$metFragCode[metFragAdducts$mode %in% 'pos']
    names(esiList) <- adductTable$adduct[metFragAdducts$mode %in% 'pos']
    } else {
    esiList <- adductTable$metFragCode[metFragAdducts$mode %in% 'neg']
    names(esiList) <- adductTable$adduct[metFragAdducts$mode %in% 'neg']
    }
    bestAnnoSub <- bestAnnoSub[bestAnnoSub$ESI_type %in% names(esiList), ]
    
    # which query spectra still have smiles
     querySpecIndxTmp <- match(bestAnnoSub$specNamesTmp, names(compSpectra(object)))  
     bestAnnoSub <- bestAnnoSub[!is.na(querySpecIndxTmp), ]
     querySpecIndxTmp <- querySpecIndxTmp[!is.na(querySpecIndxTmp)]
     featureSubSet <- unique(querySpecIndxTmp)
     querySpecDbId <- tapply(querySpecIndxTmp, paste0(bestAnnoSub$DBid, 
                                                      bestAnnoSub$ESI_type), 
                             paste0, collapse=' ')
     bestAnnoSub <- bestAnnoSub[duplicated(paste0(bestAnnoSub$DBid, bestAnnoSub$ESI_type)) == FALSE, ]
     bestAnnoSub <- cbind(bestAnnoSub, querySpecDbId=querySpecDbId[match(paste0(bestAnnoSub$DBid, bestAnnoSub$ESI_type), names(querySpecDbId))])
    
    # empty SMILES entries
    bestAnnoSub <- bestAnnoSub[bestAnnoSub$SMILES != "", , drop=FALSE]
    
    # non-covalently bound SMILES
    nonCovIndx <- grepl('\\..+', bestAnnoSub$SMILES)
    if(any(nonCovIndx)){
      message('The following ', round(sum(nonCovIndx), 0), 
           ' SMILES codes contain non-covalently bound groups and will be removed:\n', 
           paste0(which(nonCovIndx), '. ', bestAnnoSub$SMILES[nonCovIndx], 
                  collapse='\n'))
      flush.console()
      bestAnnoSub <- bestAnnoSub[nonCovIndx == FALSE, , drop=FALSE]
    }
    # writing local sdf
    message('Generating local sdf database file of ', nrow(bestAnnoSub), 
            ' unique SMILES codes...please wait\n')
    flush.console()
    
    # convert first to SMILES set object and then SDF (ChemmineR/OB)
    SDFall <- suppressWarnings(ChemmineR::smiles2sdf(bestAnnoSub$SMILES))
    SDFall@ID <- as.character(bestAnnoSub$DBid)
    validSMILES <- ChemmineR::validSDF(SDFall)
    
    if(any(validSMILES == FALSE)){
      # if more than 10 to print only print 10
      invalIndx <- ifelse(which(validSMILES == FALSE) > 10, 
                          which(validSMILES == FALSE)[1:10], 
                          which(validSMILES == FALSE))
       
      message('The following ', round(sum(validSMILES == F), 0), 
              ' SMILES code(s) produced invalid SDF files when converted using the ?smiles2sdf function:\n',
              paste0(bestAnnoSub$SMILES[invalIndx], collapse = '\n'), '...(up to first 10 shown)\n')       
      flush.console()
      SDFall <- SDFall[validSMILES]
      bestAnnoSub <- bestAnnoSub[validSMILES, , drop=FALSE]
    }
    # precursor masses
    ms1mz <- sapply(metaData(object), function(x) 
                    x[grep('_MS1_mz$', names(x))][[1]][1])
    
    if(Parameters(object)$nCores > 0){
      pmt <- proc.time()
      
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
      
      progress <- function(n){if(n %% 50 == 0){cat(paste0(n, ' of ', length(featureSubSet),
                                                          ' spectra complete.\n'))}}
      if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
      
      # foreach and dopar from foreach package
      metFragResults <- foreach(i=featureSubSet, .packages = 'ChemmineR', .options.snow=opts) %dopar% {
        dbEntryIndxTmp <- grep(paste0('^', i, ' | ', i, '$|', '^', i, '$'), 
                               bestAnnoSub$querySpecDbId)
        metFragResTmp <- metFragCl(massSpectrum=compSpectra(object)[[i]][, c('mass', 'intensity')], 
                         precMass=ms1mz[i], 
                         compSpecName=names(compSpectra(object))[i],
                         dbEntryTable=bestAnnoSub[dbEntryIndxTmp, c('WebAddress', 'DBid', 'ESI_type', 'DBname', 'SMILES'), drop=FALSE], 
                         SDFtmp=SDFall[dbEntryIndxTmp], metFragJar=metFragJar, 
                         mode=Parameters(object)$mode, 
                         frag_mzabs=frag_mzabs, 
                         esiList=esiList, maxTreeDepth=maxTreeDepth)
        # remove row.names to save on space
        row.names(metFragResTmp) <- NULL
        return(metFragResTmp)
      }
      # stop SNOW cluster
      parallel::stopCluster(cl)  
      object@inSilico$MetFrag[featureSubSet] <- metFragResults
      proc.time() - pmt
    } else { # single threaded
      if(verbose == TRUE){ pb <- txtProgressBar(max=max(featureSubSet), style=3)}
      metFragResults <- vector('list', length(featureSubSet))
    for(i in featureSubSet){
      # progress bar
      if(verbose==TRUE){setTxtProgressBar(pb, i)} 
      
      dbEntryIndxTmp <- grep(paste0('^', i, ' | ', i, '$|', '^', i, '$'), 
                             bestAnnoSub$querySpecDbId)
      metFragResTmp <- metFragCl(massSpectrum=compSpectra(object)[[i]][, c('mass', 'intensity')], 
                                       precMass=ms1mz[i], 
                                       compSpecName=names(compSpectra(object))[i],
                                       dbEntryTable=bestAnnoSub[dbEntryIndxTmp, c('WebAddress', 'DBid', 'ESI_type', 'DBname', 'SMILES'), drop=FALSE], 
                                       SDFtmp=SDFall[dbEntryIndxTmp], 
                                       metFragJar=metFragJar, 
                                       mode=Parameters(object)$mode, 
                                       frag_mzabs=frag_mzabs,
                                       esiList=esiList, maxTreeDepth=maxTreeDepth)
      # remove row.names to save on space
      row.names(metFragResTmp) <- NULL
      object@inSilico$MetFrag[[i]] <- metFragResTmp
    }  
   } # if parallel
    
    # if necessary then add to comments table
    if(autoPossId ==  TRUE){
      # message('CAUTION: do you wish to automatically add the most probable database annotations to the comments table (column possible_identity)?\n type (Y/N) and press [enter] to continue (NB. matchSpectralDB matches will not be overwritten):')
      # flush.console()
      # overWritePossId <- readline()
      # if(overWritePossId == 'Y'){
      bestCandidates <-  sapply(MetFrag(object), function(x){
          bestCand <- matrix('', ncol=2)
          if(!is.null(x)){
            if(ncol(x) > 1){
              propExTmp <- as.numeric(as.character(x[, 'propIntEx']))
              propExTmp[is.na(propExTmp)] <- 0
              scoreTmp <- x[, 'Score']
              scoreTmp[is.na(scoreTmp)] <- 0
           meanScore <- {propExTmp + scoreTmp}/2
           if(any(meanScore >= minMetFragScore)){
           bestCand <- x[which.max(meanScore), 'DBname']
           esiTypeTmp <- x[which.max(meanScore), 'ESI_type']
           bestCand <- cbind(bestCand, esiTypeTmp)
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
            message('The following automatic possible annotations have been identified as possible contaminants (i.e. occuring more than ', possContam, ' times see possContam argument) and will be named "possible_contaminant" in the Comments table and flagged "metID.metFrag":\n', paste0(1:sum(indxTmp), '. ', names(possContaminants)[indxTmp], ' (occurs n=', possContaminants[indxTmp], ' times)\n'))
            flush.console()
            contamIndx <- (bestCandidates[1, ] %in% names(possContaminants[indxTmp]))
            contamSpec <- colnames(bestCandidates)[contamIndx]
            # add to comments
            indxTmp <- metIDcomments$compSpectrum %in% contamSpec
            metIDcomments$user_comments[indxTmp] <- 'metID.metFrag'
            metIDcomments$possible_identity[indxTmp] <- 'possible_contaminant'
            bestCandidates <- bestCandidates[, contamIndx == FALSE, drop=FALSE]
          }
          
          message(ncol(bestCandidates), ' composite spectra out of ', nrow(metIDcomments), ' total (', round((ncol(bestCandidates)/nrow(metIDcomments)) * 100, 2), '%) have been annotated based on a mean metFrag score (>= ', round(minMetFragScore, 2), ').\n\n')
          flush.console()
          message('These annotations will now be added to the "metID comments" table in compMS2Explorer and flagged as "metID.metFrag".\n')
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
          metIDcomments$user_comments[!is.na(indxTmp)] <- 'metID.metFrag'
          Comments(object) <- metIDcomments
          } else {
          message('no candidates were automatically annotated by metID.metFrag try lowering the mean score.\n')
          flush.console()
          }
      }
     }
    return(object)
}) # end function
