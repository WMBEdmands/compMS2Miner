#' Quantitative Structure-Retention Relationship modelling (QSRR) using molecular descriptors and randomForest modelling
#' 
#' @details Based on the method described in Cao \emph{et. al.} \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419193/} and use the \code{\link{caret}} package (see tutorial: \url{http://topepo.github.io/caret/rfe.html} for the recursive feature selection. randomForest method utilized).
#' calculates a quantitative structure-retention relationship model
#' the default is to use the putative annotations included in the "metID comments" table of \code{\link{compMS2Explorer}} the putative annotations in the possible_identity column of the metID comments interactive table must match perfectly the database entry names found in the "best annotations" table (e.g. ensure correct matching by copy and pasting the possible compound identity in to the possible_identity column of the "metID comments" table). The metID.rtPred functions calculates molecular descriptors for all database entries in the "Best annotations" panel using the rcdk package. 
#' 
#' The molecular descriptors are then cleaned in the following sequence:
#' \enumerate{
#'  \item removing any molecular descriptors with greater than 10\% missing values.
#'  \item removing any molecular descriptors with near zero variance using the function \code{\link{nearZeroVar}} from the caret package.
#'  \item a correlation matrix of remaining molecular descriptors is calculated and
#' molecular descriptors with a standard deviation are removed.
#' \item finally any molecular descriptors with a high pair-wise correlation (>= 0.9 pearson product moment) are identified and the molecular descriptors with the largest mean absolute correlation of each group are removed. see function \code{\link{findCorrelation}} from the caret package.
#' }
#' 
#' The calculation of molecular descriptors for a large number of database entries is a potentially time-consuming process and is therefore only needs to be conducted once and the results of the process saved in the \linkS4class{CompMS2} object.
#' 
#' The caret package function \code{\link{rfe}} function is then used to identify the optimum set of remaining molecular descriptors to predict retention time. A plot should appear showing the correlation between the actual and predicted retention times of the training set.  
#' 
#' A possible workflow sequence would consist of initial examination of the results in \code{\link{compMS2Explorer}} with putative annotation of metabolites followed by use of the \code{\link{metID.rtPred}} function. After the first time the \code{\link{metID.rtPred}} function has run a new plot will appear in the \code{\link{compMS2Explorer}} gui where the "Best Annotations" closest
#' to the randomForest model predicted retention times can be easily visualized.
#' After more identifications have been made and additional putative annotations 
#' have been included in the "metID comments" table the \code{\link{metID.rtPred}}
#' function can be ran a second time. It should be much faster than the first 
#' as molecular descriptors have already been calculated and cleaned for all entries.
#' 
#' @param object A "compMS2" class object.  
#' @param standardsTable data.frame of standard compounds. The standard compounds should have been acquired using the same chromatographic method as the metabolomic dataset. If this argument is supplied then this table will be used to calculate the \code{\link{randomForest}} retention time prediction model rather than the possible_identity annotations from the "met ID comments" table. The table which must contain at mimimum the following 3 column names and an error will be returned if this is not the case (will ignore case e.g. both the column names SMILES or smiles are acceptable): 
#' \enumerate{
#' \item compound "character" type of compound names.
#' \item smiles "character" type of SMILES codes.
#' \item RT "numeric" type of retention time values (in seconds)
#' }
#' N.B. The data.frame may also contain additional columns 
#' @param descriptors character vector of molecular descriptor class names from 
#' \code{\link{get.desc.names}}. If NULL then all molecular descriptors will be considered.
#' @param removeOut logical (default = TRUE). If true outliers identified by
#' Tukey's method that is a retention time deviation of any of the training
#' set compounds greater than 1.5 * the interquartile range will be removed and
#' the QSRR model will be recalculated.
#' @param propMissing numeric maximum proportion of missing values to include a
#' molecular descriptor (values 0-1, default=0.1 i.e. maximum 10\% missing values).
#' @param propZero numeric maximum proportion of zero values to include a molecular
#' descriptor (values 0-1, default=0.2 i.e. maximum 20\% zero values).
#' @param corrPairWise numeric minimum pair-wise Pearson Product moment 
#' correlation value (values 0-1, default = 0.9),  if any molecular descriptors 
#' have high pair-wise correlation then the variables with the largest mean 
#' absolute correlation of each group are removed.
#' @param ... additional arguments to \code{\link{nearZeroVar}}.'
#' @param verbose logical if TRUE display progress bars.
#' @seealso \code{\link{nearZeroVar}}, \code{\link{rfe}}, \code{\link{randomForest}}.
#' @source Predicting retention time in hydrophilic interaction liquid chromatography mass spectrometry and its use for peak annotation in metabolomics \emph{et. al.} Metabolomics 2015 \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419193/}  
#' @examples 
#' compMS2Example <- metID(compMS2Example, 'rtPred')
#' @export
setGeneric("metID.rtPred", function(object, ...) standardGeneric("metID.rtPred"))

setMethod("metID.rtPred", signature = "compMS2", function(object, standardsTable=NULL, 
descriptors=NULL, removeOut=TRUE, propMissing=0.1, propZero=0.2, corrPairWise=0.9,
verbose=TRUE, ...){
  # error handling
  if(!require(rcdk)){
    stop('package rcdk must be installed')
  }
  if(!require(randomForest)){
    stop('package randomForest must be installed')
  }
  if(!require(caret)){
    stop('package caret must be installed')
  }
  
  stopifnot(!is.null(object))
  if(class(object) != "compMS2"){
    stop('argument object must be a "compMS2" class object')
  }
  if(length(DBanno(object)) == 0){
    stop('The function metID.dbAnnotate and metID.dbProb must be run before chemical similarity calculation')
  } 
  if(length(BestAnno(object)) == 0){
    stop('The function metID.dbProb must be run before chemical similarity calculation')  
  }

  # index which bestAnnos have entries
  emptyIndx <- sapply(BestAnno(object), nrow) > 0
  if(all(emptyIndx == FALSE)){
    stop('none of the composite spectra have any dbAnnotate matches')
  }
  # match bestAnno names to vertices
  # extract smiles codes 
  smilesDf <- do.call(rbind, BestAnno(object)[emptyIndx])
  nrowBestAnno <- sapply(BestAnno(object)[emptyIndx], nrow)
  smilesDf$specNamesTmp <- do.call(c, mapply(rep, names(BestAnno(object)[emptyIndx]),
                                             each=nrowBestAnno))
  
  # remove previous training set column if necessary
  # smilesDf$trainingSet <- NULL
  # extract retention times for each composite spectrum MS1 feature
  ms1Rt <- sapply(metaData(object)[emptyIndx],
                                    function(x) x[grep('_MS1_RT', 
                                                       names(x))][[1]][1])
  # match by comp spec id
  indxTmp <- match(smilesDf$specNamesTmp, names(ms1Rt))
  # add ms1 retention time to smiles dataframe
  smilesDf$ms1Rt <- ms1Rt[indxTmp]
  # remove substructures types will be expanded in later versions
  smilesDf <- smilesDf[smilesDf$SubStr_type == '', , drop=FALSE]
  # have molecular descriptors already been calculated?
  MD_alreadyCalc <- any(grepl('MD_', colnames(smilesDf)))
  
  if(MD_alreadyCalc == FALSE){
  # add precalculated molecular descriptors if available
  availDbs <- c('HMDB', 'LMSD', 'drugBank', 'T3DB', 'ReSpect')
  allDbIds <- do.call(c, lapply(availDbs, function(x){
    eval(parse(text=paste0('data(', x, ')')))  
    dbIdTmp <- eval(parse(text=paste0(x, '$Unique_DB_ID')))
    nrowDBtmp <- eval(parse(text=paste0('nrow(', x, ')')))
    names(dbIdTmp) <- paste0(x, '_', 1:nrowDBtmp)
    return(dbIdTmp)}))
  mdColNames <- unique(do.call(c, lapply(availDbs, function(x){
    eval(parse(text=paste0("colnames(", x, ")[grepl('^MD_', colnames(", x, "))]")))})))
  dbIdMatch <- match(smilesDf$DBid, allDbIds)
  # match indx
  if(any(!is.na(dbIdMatch))){
    # subset all ids
    allDbIds <- allDbIds[dbIdMatch[!is.na(dbIdMatch)]]
    # subset available databases
    availDbs <- availDbs[availDbs %in% unique(gsub('_.+', '', names(allDbIds)))]
    # extract MD columns
    allDbMds <- do.call(rbind, lapply(availDbs, function(x){
      indxTmp <- mdColNames %in% eval(parse(text=paste0("colnames(", x, ")"))) 
      nDBentTmp <- allDbIds[grep(paste0(x, '_'), names(allDbIds))]
      emptyMDs <- data.frame(matrix(NA, nrow=length(nDBentTmp), 
                                    ncol=length(mdColNames) + 1), stringsAsFactors = FALSE)
      colnames(emptyMDs) <- c('Unique_DB_ID', mdColNames)
      if(any(indxTmp)){
        rowNumTmp <- as.numeric(gsub(paste0(x, '_'), '', names(nDBentTmp)))
        eval(parse(text=paste0('data(', x,')'))) 
        dbTmp <- eval(parse(text=paste0('as.data.frame(', x,', stringsAsFactors=FALSE)'))) 
        emptyMDs[,  c('Unique_DB_ID', mdColNames[indxTmp])] <-  dbTmp[rowNumTmp, c('Unique_DB_ID', mdColNames[indxTmp]), drop=FALSE]
      }
      return(emptyMDs)
    }))
    # match indx
    dbIdMatch <- match(smilesDf$DBid, allDbMds$Unique_DB_ID)
    # add empty columns
    smilesDf[, mdColNames] <- 0
    smilesDf[!is.na(dbIdMatch), mdColNames] <- allDbMds[dbIdMatch[!is.na(dbIdMatch)], mdColNames]
    smilesDf$preCalcMD <- !is.na(dbIdMatch)
  } else {
    smilesDf$preCalcMD <- rep(all(mdColNames %in% colnames(smilesDf)), nrow(smilesDf)) 
  }
  }
  
 if(is.null(standardsTable)){
  # match comments table possible_identity to best anno table
  if(nrow(Comments(object)) == 0){
    stop('No comments have yet been made and therefore no possible_identity have been entered in the column.\n\nEither run the automatic chemical similarity metabolite identification (see ?metID.chemSim) or visualize with ?compMS2Explorer and add putative annotation names to the possible_identity column in the "met ID comments" tabbed panel') 
  }
  # extract comments table
  metIDcomments <- Comments(object)
  # match possible_identity to best anno table
  indxTmp <- match(smilesDf$specNamesTmp, metIDcomments$compSpectrum)
  smilesDf$possible_identity <- metIDcomments$possible_identity[indxTmp]
  # which dbNames and possible_identity names match
  smilesDf$nameMatch <- smilesDf$DBname == smilesDf$possible_identity
  if(all(smilesDf$nameMatch == FALSE) & all(smilesDf$possible_identity == '')){
    stop('No entries in the possible_identity column of the "Met ID comments" table.\n\nEither run the automatic chemical similarity metabolite identification (see ?metID.chemSim) or visualize with ?compMS2Explorer and add putative annotation names to the possible_identity column in the "met ID comments" tabbed panel')
  }
  # subset remove DBid duplicates but makesure not to remove comment table name matches 
  smilesDfSub <- smilesDf[duplicated(smilesDf$DBid) == FALSE | smilesDf$nameMatch, , drop=FALSE]
  
  message('Training set of ', length(unique(smilesDfSub$specNamesTmp[smilesDfSub$nameMatch])), ' unique annotations...\n')
  flush.console()
  
  if(any(smilesDfSub$preCalcMD == FALSE)){
    smilesAll <- smilesDfSub$SMILES 
    names(smilesAll) <- smilesDfSub$DBid
    names(smilesAll)[smilesDfSub$nameMatch] <- paste0('training', sprintf('%04d', 1:sum(smilesDfSub$nameMatch)), '_', smilesDfSub$DBid[smilesDfSub$nameMatch])
    smilesAll <- smilesAll[smilesDfSub$preCalcMD == FALSE]
    MD_alreadyCalc <- FALSE
  } else {
    MD_alreadyCalc <- TRUE  
  }
  } else {
    # check standards table contains the requisite columns
    colIndxTmp <- grepl('^compound$|^smiles$|^RT$', colnames(standardsTable),
                        ignore.case=TRUE)
    # error handling if any column names are missing
    if(sum(colIndxTmp) != 3){
      if(all(colIndxTmp == FALSE)){
      stop('Of the ', length(colIndxTmp), ' column names of the standardsTable data.frame none of the mimimum column names were identified (ignores case e.g. SMILES and smiles are both acceptable):\n\n1. compound\n2. smiles\n3. RT\n\nPlease check and try again\n')
      } else {
        stop('Of the ', length(colIndxTmp), ' column names of the standardsTable data.frame only the following mimimum column names were identified:\n\n', paste0(1:sum(colIndxTmp), '. ', colnames(standardsTable)[colIndxTmp], collapse = '\n'), '\n\nThe minimum column names must be named (ignores case e.g. SMILES and smiles are both acceptable):\n\n1. compound\n2. smiles\n3. RT\n\nPlease check and try again\n') 
      }
    }
    # remove duplicate smiles
    smilesDfSub <- smilesDf[duplicated(smilesDf$DBid) == FALSE, , drop=FALSE]
    if(MD_alreadyCalc){
      smilesAll <- standardsTable[, grep('^smiles$', colnames(standardsTable), ignore.case=TRUE)]
      names(smilesAll) <- paste0('training', sprintf('%04d', 1:nrow(standardsTable)))
    } else {
      smilesAll <- c(standardsTable[, grep('^smiles$', colnames(standardsTable), ignore.case=TRUE)], smilesDfSub$SMILES[smilesDfSub$preCalcMD == FALSE]) 
      names(smilesAll) <- c(paste0('training', sprintf('%04d', 1:nrow(standardsTable))), smilesDfSub$DBid[smilesDfSub$preCalcMD == FALSE])
    }  
   }
  message('Retention time will be predicted for ',
          nrow(smilesDfSub), ' unique SMILES codes...\n')
  flush.console()
  
  # if necessary calculate molecular descriptors this only needs to be performed once.
  if(MD_alreadyCalc == FALSE | !is.null(standardsTable)){
  # all molecular descriptor names
  if(is.null(descriptors)){
  descriptors <- rcdk::get.desc.names(type="all")
  }
  message('Evaluating molecular descriptors (n=', length(descriptors), 
          ') with rcdk package for ', length(smilesAll), ' SMILES...this is only performed once but can be time consuming.\n')
  flush.console()
  
  if(Parameters(object)$nCores > 0){
    # parallel
    if(!require(foreach)){
      stop('the foreach package must be installed to run a parallel process.\n')
    }
    # colnames 
    colNamesTmp <- colnames(HMDB)
    colNamesTmp <- colNamesTmp[grepl('^MD_', colNamesTmp)]
    colNamesTmp <- gsub('^MD_', '', colNamesTmp)
    
    pmt <- proc.time()
    nCores <- Parameters(object)$nCores
    actualCores <- parallel::detectCores()
    nCores <- ifelse(nCores >= actualCores, actualCores, nCores)
    message(paste0("Starting SNOW cluster with ", nCores, " local sockets..."))
    flush.console()
    cl <- parallel::makeCluster(nCores, outfile='')
    doSNOW::registerDoSNOW(cl)
    
    progSeq <- round({length(smilesAll) * seq(0, 1, 0.05)}, 0)
    progSeq[1] <- 1
    cat(paste0('Progress (', length(smilesAll), ' structures):\n'))
    progress <- function(n){if(n %in% progSeq){cat(paste0(round({n/length(smilesAll)} * 100, 0), '%  '))}}
    if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
    
    # foreach and dopar from foreach package
    allDescs <- foreach(smTmp=smilesAll, .packages=c('rcdk'), .combine = 'rbind', .options.snow=opts) %dopar% {
      mol <- rcdk::parse.smiles(smTmp, kekulise = FALSE)
      allDescsTmp <- tryCatch({
        suppressWarnings(eval.desc(mol, descriptors, verbose=FALSE))
      }, error=function(cond) {
        return(0)
      }, warning=function(cond){
        return(0)})
      
      if(!is.data.frame(allDescsTmp)){
        allDescsTmp <- data.frame(matrix(NA, nrow=1, ncol=length(colNamesTmp)),
                                  stringsAsFactors = FALSE)  
        colnames(allDescsTmp) <- colNamesTmp
      }
      return(allDescsTmp)
    }
    # stop SNOW cluster
    parallel::stopCluster(cl)
    proc.time() - pmt
  } else {
    pmt <- proc.time()
    if(verbose == TRUE){ pb <- txtProgressBar(max=length(smilesAll), style=3)}
    allDescs <- NULL
    for(i in 1:length(smilesAll)){
      if(verbose==TRUE){setTxtProgressBar(pb, i)}
      mol <- rcdk::parse.smiles(smilesAll[i], kekulise = FALSE)
      allDescsTmp <- tryCatch({
        suppressWarnings(eval.desc(mol, descriptors, verbose=FALSE))
      }, error=function(cond) {
        return(0)
      }, warning=function(cond){
        return(0)})
      
      if(is.data.frame(allDescsTmp)){
        if(is.null(allDescs)){
          allDescs <- data.frame(matrix(NA, nrow=length(smilesAll), ncol=ncol(allDescsTmp)),
                                 stringsAsFactors = FALSE)  
          colnames(allDescs) <- colnames(allDescsTmp)
        }
        allDescs[i, ] <- allDescsTmp
      }
    }
  }
  
  proc.time() - pmt # 2133.24 seconds 35 mins
  colnames(allDescs) <- paste0('MD_', colnames(allDescs))
  mdColNames <- colnames(allDescs)
  if(!is.null(standardsTable)){
  stdTmp <- data.frame(matrix('standards', nrow=nrow(standardsTable), 
                              ncol=ncol(smilesDfSub)), stringsAsFactors = FALSE)
  colnames(stdTmp) <- colnames(smilesDfSub)
  smilesDfSub <- rbind(stdTmp, smilesDfSub)
  if(is.null(smilesDfSub$preCalcMD)){
  smilesDfSub$preCalcMD <- TRUE
  }
  stdIndx <- grepl('standards', smilesDfSub[, 1])
  smilesDfSub$preCalcMD[stdIndx] <- FALSE
  # add std retention times SMILES and names
  smilesDfSub$DBname[stdIndx] <- standardsTable[, grepl('^compound$', colnames(standardsTable), ignore.case = TRUE)]
  smilesDfSub$SMILES[stdIndx] <- standardsTable[, grepl('^SMILES$', colnames(standardsTable), ignore.case = TRUE)]
  smilesDfSub$ms1Rt[stdIndx] <- standardsTable[, grepl('^RT$', colnames(standardsTable), ignore.case = TRUE)]
  mdColNames <- mdColNames[mdColNames %in% colnames(smilesDfSub)]
  }
  smilesDfSub[smilesDfSub$preCalcMD == FALSE, mdColNames] <- allDescs[, mdColNames]
  }
  
  indxTmp <- match(smilesDf$DBid, smilesDfSub$DBid)
  # add Moldesc to smiles
  smilesDf[!is.na(indxTmp), mdColNames] <- smilesDfSub[indxTmp[!is.na(indxTmp)], mdColNames]
  message('Cleaning molecular descriptor data (n=', length(mdColNames), ')...\n')
  flush.console()
  smilesDfColNames <- setdiff(colnames(smilesDfSub), mdColNames)
  # remove any molDesc that have too many missing values 
  indxTmp <- apply(!is.na(smilesDfSub[, mdColNames]), 2 , sum) >= (nrow(smilesDfSub) * propMissing) 
  
  message(sum(indxTmp == FALSE), ' (', sum(indxTmp), ' remaining) molecular descriptors had greater than ', round(propMissing * 100, 1), '% missing values and were removed\n ')
  flush.console()
  mdColNames <- mdColNames[indxTmp]
  
  # inpute zero to nas
  smilesDfSub[is.na(smilesDfSub)] <- 0
  
  # remove any molDesc that have too many zeros 
  indxTmp <- apply(smilesDfSub[, mdColNames] == 0, 2 , sum) <= (nrow(smilesDfSub) * propZero) 
  # message('The following ', sum(indxTmp == FALSE), ' molecular descriptors had greater than 20% missing values and were removed:\n ', paste0(colnames(allDescs)[indxTmp == FALSE], collapse = '\n'))
  message(sum(indxTmp == FALSE), ' (', sum(indxTmp), ' remaining) molecular descriptors had greater than ', round(propZero * 100, 1), '% zeros and were removed\n ')
  flush.console()
  mdColNames <- mdColNames[indxTmp]
 
  # remove near-zero variance variables
  nzv <- caret::nearZeroVar(smilesDfSub[, mdColNames], ...)
  if(length(nzv) > 0){
  mdColNames <- mdColNames[-nzv]  
  message(length(nzv), ' (', length(mdColNames), ' remaining) molecular descriptors had near-zero variance and were removed\n ')
  flush.console()
  }
  message('Calculating correlation matrix...\n')
  flush.console()
  # find highly correlation molDescs
  cor.m <- suppressWarnings(cor(apply(smilesDfSub[, mdColNames], 2, as.numeric)))
  # find rows of NAs
  noCorIndx <- !is.na(cor.m[, 1])
  message(sum(noCorIndx == FALSE), ' (', sum(noCorIndx), ' remaining) molecular descriptors had a standard deviation of zero following correlation matrix calculation and were removed.\n')
  flush.console()
  cor.m <- cor.m[noCorIndx, noCorIndx]
  mdColNames <- mdColNames[noCorIndx]
  # detect significant pair-wise comparisons
  signifPairWise <- caret::findCorrelation(cor.m, cutoff=corrPairWise)
  if(length(signifPairWise) > 0){
    mdColNames <- mdColNames[-signifPairWise]
    message(length(signifPairWise), ' (', length(mdColNames), ' remaining) molecular descriptors had a high pair-wise correlation (>= ', corrPairWise, ' pearson product moment) and the variables with the largest mean absolute correlation of each group were removed.\n')
    flush.console()  
  }
  # finally find any linear dependencies
  # comboInfo <- caret::findLinearCombos(allDescs)
    
  smilesDfSub <- smilesDfSub[, c(smilesDfColNames, mdColNames)]
 # reset trainingSet logical
  smilesDfSub$trainingSet <- FALSE
  if(!is.null(standardsTable)){
  smilesDfSub$trainingSet <- grepl('standards', smilesDfSub[, 1])  
  } else {
  indxTmp <- which(smilesDfSub$nameMatch)
  indxTmp <- indxTmp[duplicated(smilesDfSub$specNamesTmp[indxTmp]) == FALSE]
  smilesDfSub$trainingSet[indxTmp] <- TRUE
  }
  # add smiles sub to smiles df
  indxTmp <- match(smilesDf$DBid, smilesDfSub$DBid)
  addColNames <- setdiff(colnames(smilesDfSub), colnames(smilesDf))
  if(length(addColNames) > 0){
  smilesDf <- cbind(smilesDf[!is.na(indxTmp), ], smilesDfSub[indxTmp[!is.na(indxTmp)], addColNames, drop=FALSE]) 
  }
  # add training set info
  smilesDf$trainingSet <- FALSE
  indxTmp <- match(paste0(smilesDf$DBid, smilesDf$specNamesTmp), paste0(smilesDfSub$DBid, smilesDfSub$specNamesTmp))
  smilesDf$trainingSet[!is.na(indxTmp)] <- smilesDfSub$trainingSet[indxTmp[!is.na(indxTmp)]]
  ## end molecular descriptor extraction
  
  # 2. random forest based optimum variable combination selection
  ##############################################################################
  ########################### http://topepo.github.io/caret/rfe.html ###########
  ##############################################################################
  # colnames(allDescs) %in% colnames(exampleData)]
  message('performing randomForest based recursive feature elimination to identify optimal set of molecular descriptors for retention time prediction...Please wait\n')
  flush.console()
  
   rfRFE <- list(summary = defaultSummary,
                 fit = function(x, y, first, last, ...){
                   library(randomForest)
                   randomForest(x, y, importance = first, ...)
                 },
                 pred = function(object, x)  predict(object, x),
                 rank = function(object, x, y){
                   vimp <- varImp(object)
                   vimp <- vimp[order(vimp$Overall,decreasing = TRUE), , drop = FALSE]
                   vimp$var <- rownames(vimp)
                   vimp
                 },
                 selectSize = pickSizeBest,
                 selectVar = pickVars)
  ctrl <- caret::rfeControl(functions = rfRFE,
                            method = "repeatedcv",
                            repeats = 5,
                            verbose = FALSE,
                            allowParallel=ifelse(Parameters(object)$nCores > 0, TRUE, FALSE))
  
  ctrl$returnResamp <- "all"
  mdIndxTmp <- grep('MD_', colnames(smilesDfSub))
  subsets <- c(1:5, seq(10, length(mdIndxTmp), 5))
  trainingIndx <- as.logical(smilesDfSub$trainingSet)
  if(Parameters(object)$nCores > 0){
    if(!require(foreach)){
      stop('package foreach must be installed to use this function in parallel')
    }
    if(!require(doSNOW)){
      stop('package doSNOW must be installed to use this function in parallel')
    }
    message(paste0("Starting SNOW cluster with ", Parameters(object)$nCores, " local sockets..."))
    flush.console()
    cl <- parallel::makeCluster(Parameters(object)$nCores) 
    doSNOW::registerDoSNOW(cl)
   }
  trainingSet <- smilesDfSub[trainingIndx, mdIndxTmp]
  trainingSet <- apply(trainingSet, 2, as.numeric)
  rfProfile <- caret::rfe(trainingSet,
                          as.numeric(smilesDfSub$ms1Rt[trainingIndx]), sizes = subsets, 
                          rfeControl = ctrl)
  
  if(Parameters(object)$nCores > 0){
  # stop SNOW cluster
  parallel::stopCluster(cl) 
  }
  message('recursive feature selection result:\n')
  flush.console()
  print(rfProfile)
  
  # Parameters(object)$optimMolDesc <- paste0(rfProfile$optVariables, collapse='; ')
  ##############################################################################
  ########################### http://topepo.github.io/caret/rfe.html ###########
  ##############################################################################
   newDataTmp <- smilesDf[, rfProfile$optVariables, drop=FALSE]
  # inpute zero any remaining missing values
  newDataTmp[is.na(newDataTmp)] <- 0
  smilesDf$predRts <- predict(rfProfile$fit, newdata = newDataTmp)
  if(is.null(standardsTable)){
  corrCoeffRealPred <- suppressWarnings(cor(smilesDf$ms1Rt[as.logical(smilesDf$trainingSet)], smilesDf$predRts[as.logical(smilesDf$trainingSet)]))
} else {
  predRts <- predict(rfProfile$fit, newdata = trainingSet)
  corrCoeffRealPred <- cor(as.numeric(smilesDfSub$ms1Rt[as.logical(smilesDfSub$trainingSet)]), 
                           predRts)
  standardsTable$predRts <- predRts
}
  
  if(corrCoeffRealPred < 0.5){
  message('Weak correlation (', round(corrCoeffRealPred, 3), ') between actual and randomForest predicted retention times.\nConsider re-evaluation of the putative annotations in the "metID comments table" in the compMS2Explorer application.\n')
  flush.console()
  } else {
  message('Correlation coefficient between actual and randomForest predicted retention times = ', round(corrCoeffRealPred, 3), ' (Pearson product moment).\n')
  flush.console()
  }
  Parameters(object)$corrCoeffRealPred <- corrCoeffRealPred
  # plot predicted and real retention times
  if(is.null(standardsTable)){
  plot(smilesDf$ms1Rt[as.logical(smilesDf$trainingSet)][order(smilesDf$ms1Rt[as.logical(smilesDf$trainingSet)])], col='blue', xlab='elution order', ylab='retentionTime')
  points(smilesDf$predRts[as.logical(smilesDf$trainingSet)][order(smilesDf$ms1Rt[as.logical(smilesDf$trainingSet)])], col='red')
  legend('topleft', c('experimental RT (training set)', 'RandForest predicted RT (training set)'), pch=c(1, 1), pt.cex = 1, col=c('blue', 'red'))
  legend('bottomright', c('corr. coeff.: ', round(corrCoeffRealPred, 3)))
  } else {
    orderTmp <- order(as.numeric(smilesDfSub$ms1Rt[as.logical(smilesDfSub$trainingSet)]))
    plot(as.numeric(smilesDfSub$ms1Rt[as.logical(smilesDfSub$trainingSet)])[orderTmp], col='blue', xlab='elution order', ylab='retentionTime')
    points(predRts[orderTmp], col='red')
    legend('topleft', c('experimental RT (training set)', 'RandForest predicted RT (training set)'), pch=c(1, 1), pt.cex = 1, col=c('blue', 'red'))
    legend('bottomright', c('corr. coeff.: ', round(corrCoeffRealPred, 3))) 
  }
  # calculate retention time deviation from expected
  smilesDf$predRtDev <- smilesDf$predRts - smilesDf$ms1Rt
  
  # # outlier identification and removal
  # if(removeOut == TRUE){
  # # outlier detection
  # 
  # diffPredReal <- standardsTable$predRts - standardsTable$RT
  # trainNames <- standardsTable$compound
  # }
  # names(diffPredReal) <- 1:length(diffPredReal)
  # outliers <- boxplot.stats(diffPredReal)$out
  #   if(length(outliers) > 0){
  #   message(length(outliers), ' outliers were detected (by Tukeys method i.e. beyond 1.5 * IQR):\n', 
  #           paste0(names(outliers), '. ', trainNames[as.numeric(names(outliers))], ' (rtDev = ', round(outliers, 2), ').\n'), 'These will be removed and the QSRR model recalculated.\n')  
  #     flush.console()
  #     
  #   }
  # }
  
  message('Retaining final randomForest fit with optimum set of molecular descriptors and deviation of each annotation from their predicted retention time...Please wait\n')
  flush.console()
  
  # add predicted retention times back to bestAnno table
  DBanno(object)[emptyIndx] <- lapply(which(emptyIndx), function(x){
    bestAnnoTmp <- DBanno(object)[[x]]
    indxTmp <- smilesDf$specNamesTmp %in% names(DBanno(object))[x]
    smilesDfBestSub <- smilesDf[indxTmp, , drop=FALSE]
    #add tag for molecular descriptor minus DBid 
    bestAnnoTmp[, c(mdColNames, 'predRts', 'predRtDev', 'ms1Rt', 'trainingSet')] <- ''
    indxTmp <- match(bestAnnoTmp$DBid, smilesDfBestSub$DBid)
    if(any(!is.na(indxTmp))){
      bestAnnoTmp[!is.na(indxTmp), c(mdColNames, 'predRts', 'predRtDev', 'ms1Rt', 'trainingSet')] <- smilesDfBestSub[indxTmp[!is.na(indxTmp)], c(mdColNames, 'predRts', 'predRtDev', 'ms1Rt', 'trainingSet')] 
      # reorder by rt deviation
      bestAnnoTmp <- bestAnnoTmp[order(abs(as.numeric(bestAnnoTmp$predRtDev))), , drop=FALSE]
    }
    return(bestAnnoTmp)})
  # add model to object
  rtPred(object) <- list(rfModel=rfProfile, standardsTable=standardsTable)
  return(object)
}) # end function
