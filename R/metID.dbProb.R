#' selects best annotations based on substructure annotations identified
#' by \code{\link{subStructure.Annotate}} 
#' @description Most probable database annotations either automatically 
#' decided based on substructure type detected by the 
#' \code{\link{subStructure.Annotate}} or user supplied most probable annotations 
#' one composite spectrum at a time. Additionally any substructure either 
#' neutral loss or product ion with an available SMILES code will be matched
#' against all of the available annotations. The function utilizes the 
#' \code{\link{cmp.similarity}} function set to mode 2 (using the size of the 
#' descriptor intersection over the size of the smaller descriptor, to deal
#' with compounds that vary alot in size) of the ChemmineR package to calculate
#' the similarity of the substructure to the annotated structure. The average 
#' score between 0-1 of all of the substructures annotated by the 
#' \code{\link{subStructure.Annotate}}
#' function is returned in a new column in the
#' "Best Annotations" panel in the \code{\link{compMS2Explorer}} and the annotations ranked accordingly. Additionally, if 
#' either a database annotation corresponds to a substructure type annotated 
#' (e.g. glucuronide) or is the name of the database entry contains the substructure
#' name (case-insensitive) then this will be give a maximum top score of 1. 
#'  
#' 
#' @param object a compMS2 class object
#' @param nameFeat character of a unique name of a single composite spectra of
#' interest. If not supplied (default) all most probable annotations are decided automatically
#' and for all composite spectra. Previous most probable annotations will not
#' be overwritten if the function is run more than once.
#' @param DBids unique database identifier for a specific composite spectrum, in
#' combination with nameFeat argument.
#' @param minTimesId numeric (default = 2) the minimum number of times a particular
#' substructure type must be identified for it to be considered. This helps to
#' limit consideration of neutral losses/fragments that have been identified
#' once for example by chance. 
#' @param verbose logical if TRUE display progress bars. 
#' @return a compMS2 class object with most probable annotation(s) 
#' @seealso \code{\link{subStructure.Annotate}}, \code{\link{metID.dbAnnotate}},
#' \code{\link{cmp.similarity}}. 
#' @examples 
#' compMS2Example <- metID(compMS2Example, 'dbProb')
#' @export
setGeneric("metID.dbProb", function(object, ...) standardGeneric("metID.dbProb"))

setMethod("metID.dbProb", signature = "compMS2", function(object, nameFeat = NULL,
                                                          DBids = NULL, minTimesId=2,
                                                          verbose=TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if(length(DBanno(object)) == 0){
    stop("the dbAnnotate function has not yet been run the CompMS2 class object does not
         contain any potential annotations")
  } 
  if(!require(ChemmineR)){
    stop('The package ChemmineR is required to utilize this function.')
  }
  if(!require(ChemmineOB)){
    stop('The package ChemmineOB is required to utilize this function.')
  }
  
  if(!is.null(nameFeat) | !is.null(DBids)){
    if(is.null(DBids)){
      stop("A character refering to a composite spectrum to must be supplied")
    }
    if(is.null(DBids)){
      stop("A character vector containing best candidates for the feature selected 
           must be supplied")
    }
    feat.indx.tmp <- which(names(BestAnno(object)) == nameFeat)
    if(length(feat.indx.tmp) == 0){
      stop("no match was found between ", nameFeat, 
           " and any composite spectrum, please check the name and try again")
    }
    # if(length(BestAnno(object)) == 0){
    #   BestAnno.tmp <- vector("list", length(compSpectra(object)))
    #   BestAnno(object) <- BestAnno.tmp
    #   names(BestAnno(object)) <- names(compSpectra(object))
    # }
    dbAnno.tmp <- DBanno(object)[[feat.indx.tmp]]
    # subset dbannotations and add in to BestAnno object slot
    dbAnno.tmp <- dbAnno.tmp[which(dbAnno.tmp$DBid %in% DBids), , drop = FALSE]
    if(is.null(BestAnno(object)[[feat.indx.tmp]])){
      BestAnno(object)[[feat.indx.tmp]] <- dbAnno.tmp
    } else {
      dbanno.tmp <- rbind(BestAnno(object)[[feat.indx.tmp]], dbAnno.tmp)
      BestAnno(object)[[feat.indx.tmp]] <- dbAnno.tmp
    }
    
    } else {
      # automatically select the most likely candidates based on the substructures
      # identified
      # if(length(BestAnno(object)) == 0){
      #   BestAnno.tmp <- vector("list", length(compSpectra(object)))
      #   BestAnno(object) <- BestAnno.tmp
      #   names(BestAnno(object)) <- names(compSpectra(object))
      # }
      
      # indx empty db anno results
      emptyIndx <- sapply(DBanno(object), nrow) > 0
      # run in parallel if possible
      if(Parameters(object)$nCores > 0){
        # parallel
        if(!require(foreach)){
          stop('the foreach package must be installed to run a parallel process.\n')
        }
        pmt <- proc.time()
        nCores <- Parameters(object)$nCores
        actualCores <- parallel::detectCores()
        nCores <- ifelse(nCores >= actualCores, actualCores, nCores)
        message(paste0("Starting SNOW cluster with ", nCores, " local sockets..."))
        flush.console()
        cl <- parallel::makeCluster(nCores, outfile='')
        doSNOW::registerDoSNOW(cl)
        
        progress <- function(n){if(n %% 50 == 0){cat(paste0(n, ' of ', sum(emptyIndx),
                                                            ' spectra complete.\n'))}}
       
        if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}

        # foreach and dopar from foreach package
        Best.Anno.tmp <- foreach(compSpec=compSpectra(object)[emptyIndx], dbAnno.tmp=DBanno(object)[emptyIndx], .packages=c('ChemmineR', 'ChemmineOB'), .options.snow=opts) %dopar% {
          SubStr.types <- ''
          SubStr.SMI <- character()
          if(!is.null(compSpec)){
            SubStr.types <- unlist(compSpec[, c("Frag.ID.type", 
                                                "interfrag.loss.type",
                                                "Neutral.loss.type")])
            SubStr.types <- unlist(strsplit(ifelse(SubStr.types == '', 'noID', SubStr.types), ";")) 
            
            # which(as.matrix(compSpec) %in% SubStr.types[2], arr.ind = TRUE)
            SubStr.SMI <- unlist(compSpec[, c("Frag.ID.SMILES", 
                                              "interfrag.loss.SMILES",
                                              "Neutral.loss.SMILES")])
            SubStr.SMI <- unlist(strsplit(ifelse(SubStr.SMI == '', 'noID', SubStr.SMI), ";"))
            # remove any substructure id'ed more than once
            dupIndx <- duplicated(SubStr.SMI) == FALSE
            SubStr.types <- SubStr.types[dupIndx]
            SubStr.SMI <- SubStr.SMI[dupIndx]
            # minimum times identified
            # above minimum number of occurences
            nTimesIdent <- table(SubStr.types)
            nTimesIdent <- nTimesIdent >= minTimesId & names(nTimesIdent) != 'noID'
            nTimesIdent <- SubStr.types %in% names(nTimesIdent)[nTimesIdent]
            SMI.avail <- SubStr.SMI != 'noSMILES' & nTimesIdent
            SubStr.types <- c('', SubStr.types[SMI.avail])
            SubStr.SMI <- SubStr.SMI[SMI.avail]
          }
          # give all of the annotations with a substructure type identified a top score
          # of 1
          dbAnno.tmp$subStrScore <- 0
          subSetIndx <- dbAnno.tmp$SubStr_type %in% SubStr.types
          if(any(subSetIndx)){
            # subset substructures or the dbname contains the substr type name
            dbAnno.tmp$BestAnno <- subSetIndx
            subStrIndx <- dbAnno.tmp$SubStr_type != '' & subSetIndx
            subStrInName <- paste0(SubStr.types[SubStr.types != ''], collapse='|')
            if(subStrInName != ''){
              subStrIndx[grepl(subStrInName, dbAnno.tmp$DBname, ignore.case = TRUE) & subSetIndx] <- TRUE
            }
            dbAnno.tmp$subStrScore[subStrIndx] <- 1
            
            # identify presence of any annotated substructures in SMILES
            remainingIndx <- subStrIndx == FALSE & subSetIndx == TRUE
            if(length(SubStr.SMI) > 0 & any(remainingIndx)){
              sdfSubStr <- try(suppressWarnings(ChemmineR::smiles2sdf(SubStr.SMI)))
              if(class(sdfSubStr) != 'try-error'){
                apSubStr.SMI <- ChemmineR::sdf2ap(sdfSubStr)
                apDbAnno.SMI <- ChemmineR::sdf2ap(suppressWarnings(ChemmineR::smiles2sdf(dbAnno.tmp$SMILES[remainingIndx])))
                simScoreTmp <- vector('numeric', length(apDbAnno.SMI))
                for(j in 1:length(apSubStr.SMI)){
                  for(k in 1:length(apDbAnno.SMI)){
                    # if the functional group present then score 1
                    fgPres <- cmp.similarity(apSubStr.SMI[j], apDbAnno.SMI[k], mode=2)
                    simScoreTmp[k] <- simScoreTmp[k] + fgPres
                  }
                }
                # average all similarity scores
                dbAnno.tmp$subStrScore[remainingIndx] <- simScoreTmp/length(apSubStr.SMI)
                dbAnno.tmp <- dbAnno.tmp[order(dbAnno.tmp$subStrScore, decreasing=TRUE), , drop=FALSE]
              }
            }
          }
            return(dbAnno.tmp)
        }
        # stop SNOW cluster
        parallel::stopCluster(cl)
        proc.time() - pmt
      } else {
      pmt <- proc.time()
      if(verbose == TRUE){ pb <- txtProgressBar(max=length(DBanno(object)[emptyIndx]), style=3)}
      Best.Anno.tmp <- vector('list', length(DBanno(object)[emptyIndx]))
      for(i in 1:length(DBanno(object)[emptyIndx])){
        if(verbose==TRUE){setTxtProgressBar(pb, i)}
        SubStr.types <- ''
        SubStr.SMI <- character()
        compSpec <- compSpectra(object)[emptyIndx][[i]]
        dbAnno.tmp <- DBanno(object)[emptyIndx][[i]]
        if(!is.null(compSpec)){
        SubStr.types <- unlist(compSpec[, c("Frag.ID.type", 
                                                            "interfrag.loss.type",
                                                            "Neutral.loss.type")])
        SubStr.types <- unlist(strsplit(ifelse(SubStr.types == '', 'noID', SubStr.types), ";")) 
        
        # which(as.matrix(compSpec) %in% SubStr.types[2], arr.ind = TRUE)
        SubStr.SMI <- unlist(compSpec[, c("Frag.ID.SMILES", 
                                                                  "interfrag.loss.SMILES",
                                                                  "Neutral.loss.SMILES")])
        SubStr.SMI <- unlist(strsplit(ifelse(SubStr.SMI == '', 'noID', SubStr.SMI), ";"))
        # remove any substructure id'ed more than once
        dupIndx <- duplicated(SubStr.SMI) == FALSE
        SubStr.types <- SubStr.types[dupIndx]
        SubStr.SMI <- SubStr.SMI[dupIndx]
        # minimum times identified
        # above minimum number of occurences
        nTimesIdent <- table(SubStr.types)
        nTimesIdent <- nTimesIdent >= minTimesId & names(nTimesIdent) != 'noID'
        nTimesIdent <- SubStr.types %in% names(nTimesIdent)[nTimesIdent]
        SMI.avail <- SubStr.SMI != 'noSMILES' & nTimesIdent
        SubStr.types <- c('', SubStr.types[SMI.avail])
        SubStr.SMI <- SubStr.SMI[SMI.avail]
        }
        # give all of the annotations with a substructure type identified a top score
        # of 1
        dbAnno.tmp$subStrScore <- 0
        subSetIndx <- dbAnno.tmp$SubStr_type %in% SubStr.types
        
        if(any(subSetIndx)){
        # subset substructures or the dbname contains the substr type name
        dbAnno.tmp$BestAnno <- subSetIndx
        subStrIndx <- dbAnno.tmp$SubStr_type != '' & subSetIndx
        subStrInName <- paste0(SubStr.types[SubStr.types != ''], collapse='|')
        if(subStrInName != ''){
        subStrIndx[grepl(subStrInName, dbAnno.tmp$DBname, ignore.case = TRUE) & subSetIndx] <- TRUE
        }
        dbAnno.tmp$subStrScore[subStrIndx] <- 1
        
        # identify presence of any annotated substructures in SMILES
        remainingIndx <- subStrIndx == FALSE & subSetIndx == TRUE
        if(length(SubStr.SMI) > 0 & any(remainingIndx)){
          sdfSubStr <- try(suppressWarnings(ChemmineR::smiles2sdf(SubStr.SMI)))
          if(class(sdfSubStr) != 'try-error'){
           apSubStr.SMI <- ChemmineR::sdf2ap(sdfSubStr)
           apDbAnno.SMI <- ChemmineR::sdf2ap(suppressWarnings(ChemmineR::smiles2sdf(dbAnno.tmp$SMILES[remainingIndx])))
        simScoreTmp <- vector('numeric', length(apDbAnno.SMI))
        for(j in 1:length(apSubStr.SMI)){
          for(k in 1:length(apDbAnno.SMI)){
          # if the functional group present then score 1
          fgPres <- cmp.similarity(apSubStr.SMI[j], apDbAnno.SMI[k], mode=2)
          simScoreTmp[k] <- simScoreTmp[k] + fgPres
          }
        }
        # average all similarity scores
        dbAnno.tmp$subStrScore[remainingIndx] <- simScoreTmp/length(apSubStr.SMI)
        dbAnno.tmp <- dbAnno.tmp[order(dbAnno.tmp$subStrScore, decreasing=TRUE), , drop=FALSE]
          }
         }
        }
       Best.Anno.tmp[[i]] <- dbAnno.tmp
       }
      }
      proc.time() - pmt
      names.tmp <- names(compSpectra(object)[emptyIndx])
      names(Best.Anno.tmp) <- names.tmp
      DBanno(object)[names.tmp] <- Best.Anno.tmp
      # CHECK   IF ANY  are still unannotated and add to comments table
      noAnno <- sapply(BestAnno(object), nrow) == 0
      metIDcomments <- Comments(object)
      noAnnoComm <- metIDcomments$possible_identity == 'no_annotations'
      metIDcomments$possible_identity[noAnno == FALSE & noAnnoComm] <- ''
      emptyComment <- metIDcomments$possible_identity == ''  
      metIDcomments$possible_identity[noAnno & {noAnnoComm | emptyComment}] <- 'no_annotations'
      Comments(object) <- metIDcomments
      # names(BestAnno(object)) <- names.tmp
      return(object)  
}}) # end function
