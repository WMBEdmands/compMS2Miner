#' Constructor for compMS2 class object from a peak table and MS2 mzXML/mzML/.mgf
#' file(s)
#'
#' @description Matches MS1 features to MS2 spectra (.mzXML/.mzML/.mgf) files 
#' based on a mass-to-charge and retention time tolerance. Composite spectra 
#' and other data can subsequently be visualized during any stage of the compMS2Miner
#' processing workflow. Composite spectra can be denoised, ion signals grouped 
#' and summed, substructure groups identified, common Phase II metabolites
#' predicted and features matched to data bases monoisotopic mass data 
#' and insilico MS2 fragmentation data.
#' The resulting data can then be readily curated by sending to a local or online
#' couchDB database. 
#'
#' @param MS1features either a data.frame, full file path as a character string to a  .csv file of a MS1 feature table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows, the first 3 columns must consist of:
#' \enumerate{
#'  \item EIC number or unique peak identifier.
#'  \item mass-to-charge ratio of peak group.
#'  \item median/ peak apex retention time in seconds. 
#'  }
#'  
#'  optionally the 4th column of the MS1feature table may contain any adducts
#'  and isotopes identified by for example the CAMERA R package.
#'  If this column is present the adducts will be incorporated in to the compMS2
#'  class object and used to guide the subsequent \code{\link{metID.dbAnnotate}}
#'  function. This can be very useful for narrowing possible annotations in
#'  subsequent stages of the compMS2miner workflow and particularly in reduction
#'  of false positives annotations. The adduct annotations must
#'  consist of the following notation style for example [M-H]-, [2M+2CH3OH]2-, 
#'  [M-H+C2H4O2+Na]-. Abbreviations such as Hac (CH3COOH) for acetic acid 
#'  and ACN (i.e. C2H3N) for acetonitrile 
#'  are not acceptable formulae must be used to determine the correct
#'  elemental composition must be included. As is typical of the output of
#'  CAMERA for example multiple possible adducts can appear for the same feature
#'  where they have shared/similar expected masses.
#'    
#'  If argument is not supplied a GUI (tcltk) file selection window will open and a .csv file can then be selected. 
#' @param msDataDir character full path to a directory containing LC-MS/MS data files
#' in either the open framework .mzXML or newer .mzML file types also mascot generic format files (.mgf). If argument is
#' not supplied a GUI (tcltk) file selection window will open and the directory 
#' can be selected.
#' @param MS2files character vector of full paths to ms2 files (either .mzML, .mzXML or .mgf).
#' Alternative to choosing the directory. In this way particular files within a directory
#' or files from multiple directory locations can be specified.
#' @param nCores numeric Number of cores for parallel computation.
#' @param mode character Ionisation polarity must be either 'pos' or 'neg'.
#' @param precursorPpm numeric Parts per million mass accuracy to match MS1 features to MS2 spectra (ppm) 
#' @param ret numeric retention time tolerance to match MS1 features to MS2 spectra (+/- seconds). 
#' @param TICfilter numeric Minimum Total Ion Current to consider an MS2 spectrum. Any MS2 scan
#' below this threshold will not be considered. 
#' @param minPeaks. minimum number of fragment ions for a spectrum to be 
#' considered (default = 1).
#' @param isoWid numeric isolation width of DDA precursor ions, utilized to identify
#' potentially chimeric spectra.
#' @param verbose logical if TRUE display progress bars.
#' @return A compMS2 object                   
#' @export
compMS2Construct <-  function(MS1features = NULL, msDataDir = NULL, MS2files=NULL, 
                              nCores = NULL, 
                              argsCorrNetwork=list(obsNames=NULL, corrThresh=0.6,
                                                   corrMethod="spearman", 
                                                   delta=0.05, MTC="BH"), 
                              mode = "pos", precursorPpm = 10, ret = 10, 
                              TICfilter = 10000, minPeaks=1, isoWid=4, 
                              verbose=TRUE){
  
  message("creating compMS2 object in ", ifelse(mode == "pos","positive", 
                                                "negative"), 
          " ionisation mode")
  flush.console()
  if(!is.null(nCores)){
    if(!require(foreach)){
      stop('foreach package must be installed for parallel computation.\n')
    }
  }
  # new compMS2 object 
  object <- new("compMS2")
  # set global options
  options(stringsAsFactors = FALSE)
  # if is.null msDataDir select .mzXML/.mzML/.mgf file containing directory
  if(is.null(MS2files)){
  if(is.null(msDataDir)){
  message("Select your .mzXML/.mzML/.mgf data directory")
  flush.console()
  msDataDir <- tcltk::tk_choose.dir(default = "",  
                                    caption = "1. Select your .mzXML/.mzML/.mgf data directory")
  }
  
  # identify all mzXML/.mzML files in raw-data directory
  MS2files <- list.files(path=msDataDir, pattern = "\\.mzXML$|\\.mzML$|\\.mgf$", 
                         full.names=TRUE)
  }
  fileTypeTmp <- unique(gsub('.+\\.', '', basename(MS2files)))
  if(length(fileTypeTmp) > 1){
    stop('only 1 file type is permitted in the MS data directory (msDataDir):\n',
         'The following file types were found:\n', paste0('.', fileTypeTmp, '\n'))
  }
  message(length(MS2files), 
         paste0(" MS2 (.", fileTypeTmp, ") file(s) were detected within the directory..."))
  flush.console()
  adducts <- FALSE
  # if file type is mgf
  if(fileTypeTmp == 'mgf'){
    # file names
    fileNames <- gsub('\\.mgf$', '', basename(MS2files))
    # remove hypens
    fileNames <- gsub('-', '_', fileNames)
    # spectra
    compSpectra(object) <- vector('list', length(MS2files))
    names(compSpectra(object)) <- fileNames
    # metaData
    metaData(object) <- vector('list', length(MS2files))
    names(metaData(object)) <- fileNames
    message('reading mgf files...\n')
    flush.console()
    if(verbose == TRUE){ pb <- txtProgressBar(max=length(MS2files), style = 3)}
    for(i in 1:length(MS2files)){
      if(verbose==TRUE){setTxtProgressBar(pb, i)}
      mgfFile <- readLines(MS2files[i])
      # spectrum
      spectrumTmp <- mgfFile[grep('^[0-9]', mgfFile)]
      spectrumTmp <- do.call(rbind, strsplit(spectrumTmp, '\\t'))
      if(nrow(spectrumTmp) > 1){
      spectrumTmp <- data.frame(apply(spectrumTmp, 2, as.numeric), 
                                stringsAsFactors=FALSE)
      } else {
        spectrumTmp <- as.numeric(spectrumTmp[1, ])
        spectrumTmp <- data.frame(spectrumTmp[1], spectrumTmp[2], stringsAsFactors=FALSE)
      }
      colnames(spectrumTmp) <- c('mass', 'intensity')
      compSpectra(object)[[i]] <- spectrumTmp
      # metaData
      massIndx <- grep('PEPMASS=', mgfFile) 
      rtIndx <- grep('RTINSECONDS=', mgfFile) 
      metaDataTmp <- list(as.numeric(gsub('PEPMASS=', '', mgfFile[massIndx])), 
                          as.numeric(gsub('RTINSECONDS=', '', mgfFile[rtIndx])),
                          0)
      names(metaDataTmp) <- paste0(fileNames[i], c('_MS1_mz', '_MS1_RT', '_precursorIntensity'))
      metaData(object)[[i]] <- metaDataTmp
    }
  } else { # if .mzXML or .mzML
    # if is.null MS1features select MS1 feature table file
    if(is.null(MS1features)){
      message("Select your MS1 feature (.csv) file")
      flush.console()
      
      MS1features <- tcltk::tclvalue(tcltk::tkgetOpenFile(
        initialdir = msDataDir,
        filetypes = "{{Comma delimited} {.csv}} {{All files} *}", 
        title = "2. select the MS1 features (.csv) file you wish to match"))
    }
    
    if(is.character(MS1features)){
      # read in MS1features
      message("Reading MS1 feature table...")
      flush.console()
      MS1features <- as.data.frame(data.table::fread(MS1features, sep=",", 
                                                     header=TRUE, stringsAsFactors=FALSE))
      message("...Done")
      flush.console()
    }
    
    if(!is.data.frame(MS1features)){
      stop('The MS1features object is not a data.frame.')
    }
    # error handling
    if(!is.integer(MS1features[, 1])){
      stop('The first column of the MS1 feature table must be an integer (EIC/unique number).')
    }
    if(!is.numeric(MS1features[, 2])){
      stop('The second column of the MS1 feature table must be a numeric (mass-to-charge ratio).')
    }
    if(!is.numeric(MS1features[, 3])){
      stop('The third column of the MS1 feature table must be a numeric (retention time in seconds).')
    }
    # sort dataframe by unique ID/ EIC number
    MS1features <- MS1features[order(MS1features[, 1]), ]
    row.names(MS1features) <- seq(1, nrow(MS1features), 1)
    # see if the 4 th column contains adduct information
    adducts <- any(grepl('M\\+|M\\-', MS1features[, 4]))
    if(adducts == TRUE){
    message('Adducts/fragments were detected in the 4th column of the MS1features table. These may be used for subsequent stages of the workflow.\n')
    flush.console()
    }
  if(!is.null(nCores)){
    if(!require(foreach)){
      stop('package foreach must be installed to use this function in parallel')
    }
    if(!require(doSNOW)){
      stop('package doSNOW must be installed to use this function in parallel')
    }
    message(paste0("Starting SNOW cluster with ", nCores, " local sockets..."))
    flush.console()
    cl <- parallel::makeCluster(nCores, outfile='')
    doSNOW::registerDoSNOW(cl)
    
    message("matching MS1 peak table features to the following MS2 files: ")
    flush.console()
    mS2message <- sapply(basename(MS2files), function(x){ 
      message(x)
      flush.console()})
    
    progress <- function(n) cat(paste0(n, ' of ', length(MS2files),
                                       ' complete (', basename(MS2files)[n],
                                       ').\n'))
    if(verbose == TRUE){opts <- list(progress=progress)} else {opts <- list(progress=NULL)}
    # foreach and dopar from foreach package
    Results <- foreach(i=1:length(MS2files), .packages=c('mzR'),
                       .options.snow=opts) %dopar% {
      compMS2Create(MS2file=MS2files[i], MS1features=MS1features,
                    TICfilter=TICfilter,  precursorPpm=precursorPpm, 
                    ret=ret, adducts=adducts, isoWid=isoWid)
    }
    # stop SNOW cluster
    parallel::stopCluster(cl) 
  } else { 
    # create list to store results
    Results <- vector("list", length(MS2files))
    
    for(i in 1:length(MS2files)){
    Res.tmp <- compMS2Create(MS2file=MS2files[i], TICfilter=TICfilter, 
                             MS1features=MS1features, precursorPpm=precursorPpm,
                             ret=ret, adducts=adducts, isoWid=isoWid)
    Results[[i]] <- Res.tmp
    }
  }
  
  Results <- unlist(Results, recursive = FALSE)
  compSpectra(object) <- lapply(Results, function(x) x$spectra)
  metaData(object) <- lapply(Results, function(x) x$metaData)
  # check all compSpectra are not empty
  indxTmp <- sapply(compSpectra(object), nrow) > minPeaks
  if(any(indxTmp == FALSE)){
  compSpectra(object) <- compSpectra(object)[indxTmp]
  metaData(object) <- metaData(object)[indxTmp]
  }
  # MS1 feature table
  # MS1features(object) <- MS1features
  # if obsNames supplied then perform corrNetwork 
  if(!is.null(argsCorrNetwork$obsNames)){
    argsCorrNetwork[['object']] <- object
    argsCorrNetwork[['peakTable']] <- MS1features
    object <- do.call(metID.corrNetwork, argsCorrNetwork)  
  } else {
  message('\n"obsNames" in argsCorrNetwork argument missing. Not performing correlation network calculation.\n')
  flush.console()
  }
  }
  # add file paths and parameters
  filePaths(object) <- MS2files
  Parameters(object) <- data.frame(nCores=ifelse(is.null(nCores), 0, nCores),
                                   mode=mode, precursorPpm=precursorPpm,
                                   ret=ret, TICfilter=TICfilter, 
                                   fileType=fileTypeTmp, adducts=adducts, 
                                   isoWid=isoWid, stringsAsFactors=FALSE)
  return(object)
  } # end CompMS2obj function

setMethod("show", "compMS2", function(object) {
  if(length(filePaths(object)) > 0){
    cat("A \"compMS2\" class object derived from", length(filePaths(object)), 
        "MS2 files \n\n")
    # any non-ms2 matched
    noMs2Idx <- grepl('noMS2', names(metaData(object)))
    
    indxTmp <- noMs2Idx == FALSE
    
    Acc.tmp <- sapply(metaData(object)[which(indxTmp)], function(x){ 
      c(mean(abs(unlist(x[grep("ppmDiff", names(x))]))),
        mean(abs(unlist(x[grep("rtDiff", names(x))]))),
      length(unlist(x[grep("rtDiff", names(x))])))})
    ionsCount <- sum(sapply(compSpectra(object)[indxTmp], function(x) nrow(x)))
    
    spec.names <- names(compSpectra(object))[indxTmp]
    cat(length(unique(gsub(".+_", "", spec.names))), 
        "MS1 features were matched to", sum(Acc.tmp[3, ]), "MS2 precursor scans\n")
    cat("containing",ionsCount, "ion features\n\n")    
    
    cat("Average ppm match accuracy:",
        round(mean(Acc.tmp[1, ]), digits = 3) ,"\n")
    cat("with a ppm mass accuracy tolerance of (+/-)", Parameters(object)$precursorPpm,
        "\n\n")
    
    cat("Average retention time match accuracy:",
        round(mean(Acc.tmp[2, ]), digits = 2) ,"seconds\n")
    cat("with a retention time tolerance of (+/-)", Parameters(object)$ret,
        " seconds\n\n")
    # poss chimeric spectra
    nChimScans <- sum(sapply(metaData(object)[which(indxTmp)], function(x){ 
      idxTmp <- duplicated(unlist(x[grep('precursorScanNum', names(x))])) == FALSE
      sum(as.logical(unlist(x[grep('possChim', names(x))])[idxTmp]))}), na.rm=TRUE)
    allPrecScans <- sum(do.call(c, sapply(metaData(object)[which(indxTmp)], function(x){ 
      duplicated(unlist(x[grep('precursorScanNum', names(x))])) == FALSE})), na.rm=TRUE)
  
    cat(paste0(nChimScans, ' spectra of ', allPrecScans, ' precursor scans were identified as potentially chimeric (', round({nChimScans/allPrecScans} * 100, 2),'% an isolation width of ',  Parameters(object)$isoWid, ' Da was supplied to "compMS2Contruct").\n\n'))
    if(any(noMs2Idx)){
      cat(paste0('Additionally ', sum(noMs2Idx), ' EICs (', round({sum(noMs2Idx)/length(noMs2Idx)} * 100, 2),'%) are not matched to any MS2 scans and were added by the metID.corrNetwork function (correlation threshold >= ', Parameters(object)$corrThresh, '.\n'))
    }
    # maxInterInt <- do.call(c, lapply(metaData(object), function(x) x[duplicated(x[, 'precursorScanNum']) == FALSE, 'maxInterIons']))
    # rtPrecScans <- do.call(c, lapply(metaData(object), function(x) x[duplicated(x[, 'precursorScanNum']) == FALSE, 'retentionTime']))
    # maxInterInt <- as.numeric(maxInterInt)
    # rtPrecScans <- as.numeric(rtPrecScans)
    # plot(rtPrecScans, maxInterInt, type='h')
    # boxplot(log(maxInterInt))
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  } else {
    cat("A new empty\"compMS2\" class object")
  }
})


# split method e.g. cut(1:length(x@compSpectra), 2)
# setMethod("split", "compMS2", function(x, f){
# stopifnot(class(x) == 'compMS2')
#   
#   splitX <- vector('list', nlevels(f))
#   for(i in 1:nlevels(f)){
#   tmpCompMS2 <- new('compMS2')
#   indxTmp <- which(f == levels(f)[i])
#   tmpCompMS2@compSpectra <- x@compSpectra[indxTmp]
#   tmpCompMS2@metaData <- x@metaData
#   # tmpCompMS2@MS1features <- x@MS1features
#   tmpCompMS2@DBanno <- x@DBanno[indxTmp]
#   tmpCompMS2@BestAnno <- x@BestAnno[indxTmp]
#   subStrIndx <- x@subStrAnno$compSpecName %in% names(tmpCompMS2@compSpectra)[indxTmp]
#   tmpCompMS2@subStrAnno <- x@subStrAnno[subStrIndx, , drop=FALSE]
#   tmpCompMS2@Comments<- x@Comments[indxTmp]
#   tmpCompMS2@filePaths <- x@file.paths
#   tmpCompMS2@Parameters <- x@Parameters
#   tmpCompMS2@couchDBconn <- x@couchDBconn
#   tmpCompMS2@inSilico$MetFrag <- x@inSilico$MetFrag[indxTmp]
#   tmpCompMS2@inSilico$CFM <- x@inSilico$CFM[indxTmp]
#   splitX[[i]] <- tmpCompMS2
#   names(splitX)[i] <- i
#   }
#   return(splitX)
# })
#' @export
setGeneric("filePaths", function(object, ...) standardGeneric("filePaths"))

setMethod("filePaths", "compMS2", function(object) object@filePaths)
#' @export
setGeneric("filePaths<-", function(object, value) standardGeneric("filePaths<-"))

setReplaceMethod("filePaths", "compMS2", function(object, value){
  object@filePaths <- value
  return(object)
})

#' @export
setGeneric("compSpectra", function(object, ...) standardGeneric("compSpectra"))

setMethod("compSpectra", "compMS2", function(object) object@compSpectra)
#' @export
setGeneric("compSpectra<-", function(object, value) standardGeneric("compSpectra<-"))

setReplaceMethod("compSpectra", "compMS2", function(object, value){
  object@compSpectra <- value
  return(object)
})
#' @export
setGeneric("metaData", function(object, ...) standardGeneric("metaData"))

setMethod("metaData", "compMS2", function(object) object@metaData)
#' @export
setGeneric("metaData<-", function(object, value) standardGeneric("metaData<-"))

setReplaceMethod("metaData", "compMS2", function(object, value){
  object@metaData  <- value
  return(object)
})
#' @export
setGeneric("DBanno", function(object, ...) standardGeneric("DBanno"))

setMethod("DBanno", "compMS2", function(object) object@DBanno)
#' @export
setGeneric("DBanno<-", function(object, value) standardGeneric("DBanno<-"))

setReplaceMethod("DBanno", "compMS2", function(object, value) {
  object@DBanno <- value
  return(object)
})
#' @export
setGeneric("BestAnno", function(object, ...) standardGeneric("BestAnno"))

setMethod("BestAnno", "compMS2", function(object){ 
  dbAnnoTmp <- object@DBanno
  bestAnnoTmp <- lapply(dbAnnoTmp, function(x){ 
                        if('BestAnno' %in% colnames(x)){
                        return(x[x[, 'BestAnno'], , drop=FALSE]) 
                        }})
  return(bestAnnoTmp)
  })
#' @export
setGeneric("BestAnno<-", function(object, value) standardGeneric("BestAnno<-"))

setReplaceMethod("BestAnno", "compMS2", function(object, value){
  # check if dbAnno is empty
  if(length(object@DBanno) == 0){
    stop('There are no annotations in the DBanno slot run metID.dbAnnotate first.')
  }
  # check all names match
  indxTmp <- match(names(value), names(object@DBanno))
  if(any(is.na(indxTmp))){
    stop('BestAnno: The following names do not match:\n', paste0(names(value)[is.na(indxTmp)], collapse='\n'),
         '\nPlease check and try again.\n')
  }
  # add logical to db anno table
  for(i in 1:length(indxTmp)){
  object@DBanno[[indxTmp[i]]]$BestAnno <- object@DBanno[[indxTmp[i]]]$DBid %in% value[[i]]$DBid
  dbIdIndx <- match(object@DBanno[[indxTmp[i]]]$DBid, value[[i]]$DBid)
  # additional columns
  newCols <- setdiff(colnames(value[[i]]), colnames(object@DBanno[[indxTmp[i]]]))
  if(length(newCols) > 0){
   object@DBanno[[indxTmp[i]]][, newCols] <- 0
   object@DBanno[[indxTmp[i]]][!is.na(dbIdIndx), newCols] <- value[[i]][dbIdIndx[!is.na(dbIdIndx)], newCols]
  }
  object@DBanno[[indxTmp[i]]][!is.na(dbIdIndx), colnames(value[[i]])] <- value[[i]][dbIdIndx[!is.na(dbIdIndx)], ]
  }
  return(object)
})
#' @export
setGeneric("subStrAnno", function(object, ...) standardGeneric("subStrAnno"))

setMethod("subStrAnno", "compMS2", function(object) object@subStrAnno)
#' @export
setGeneric("subStrAnno<-", function(object, value) standardGeneric("subStrAnno<-"))

setReplaceMethod("subStrAnno", "compMS2", function(object, value){
  object@subStrAnno <- value
  return(object)
})
#' @export
setGeneric("MetFrag", function(object, ...) standardGeneric("MetFrag"))

setMethod("MetFrag", "compMS2", function(object) object@inSilico$MetFrag)
#' @export
setGeneric("MetFrag<-", function(object, value) standardGeneric("MetFrag<-"))

setReplaceMethod("MetFrag", "compMS2", function(object, value) {
  object@inSilico$MetFrag <- value
  return(object)
})
#' @export
setGeneric("CFM", function(object, ...) standardGeneric("CFM"))

setMethod("CFM", "compMS2", function(object) object@inSilico$CFM)
#' @export
setGeneric("CFM<-", function(object, value) standardGeneric("CFM<-"))

setReplaceMethod("CFM", "compMS2", function(object, value) {
  object@inSilico$CFM <- value
  return(object)
})
#' @export
setGeneric("Comments", function(object, ...) standardGeneric("Comments"))

setMethod("Comments", "compMS2", function(object){ 
  if(nrow(object@Comments) == 0){
  object@Comments <-  data.frame(compSpectrum=names(compSpectra(object)), 
                                 possible_identity=rep('', length(compSpectra(object))), 
                                 ESI_type=rep('', length(compSpectra(object))),
                                 compound_class=rep('', length(compSpectra(object))), 
                                 user_comments=rep('', length(compSpectra(object))), 
                                 stringsAsFactors = FALSE)  
  }
  object@Comments})
#' @export
setGeneric("Comments<-", function(object, value) standardGeneric("Comments<-"))


setReplaceMethod("Comments", "compMS2", function(object, value) {
  object@Comments <- value
  return(object)
})
#' @export
setGeneric("Parameters", function(object, ...) standardGeneric("Parameters"))

setMethod("Parameters", "compMS2", function(object) object@Parameters)
#' @export
setGeneric("Parameters<-", function(object, value) standardGeneric("Parameters<-"))

setReplaceMethod("Parameters", "compMS2", function(object, value) {
  
  object@Parameters <- value
  return(object)
})
#' @export
setGeneric("couchDBconn", function(object, ...) standardGeneric("couchDBconn"))

setMethod("couchDBconn", "compMS2", function(object) object@couchDBconn)
#' @export
setGeneric("couchDBconn<-", function(object, value) standardGeneric("couchDBconn<-"))

setReplaceMethod("couchDBconn", "compMS2", function(object, value) {
  
  object@couchDBconn <- value
  return(object)
})

#' @export
setGeneric("network", function(object, ...) standardGeneric("network"))

setMethod("network", "compMS2", function(object) object@network)
#' @export
setGeneric("network<-", function(object, value) standardGeneric("network<-"))

setReplaceMethod("network", "compMS2", function(object, value) {
  
  object@network <- value
  return(object)
})

#' @export
setGeneric("rtPred", function(object, ...) standardGeneric("rtPred"))

setMethod("rtPred", "compMS2", function(object) object@rtPred)
#' @export
setGeneric("rtPred<-", function(object, value) standardGeneric("rtPred<-"))

setReplaceMethod("rtPred", "compMS2", function(object, value) {
  
  object@rtPred <- value
  return(object)
})
#' @export
setGeneric("spectralDB", function(object, ...) standardGeneric("spectralDB"))

setMethod("spectralDB", "compMS2", function(object) object@spectralDB)
#' @export
setGeneric("spectralDB<-", function(object, value) standardGeneric("spectralDB<-"))

setReplaceMethod("spectralDB", "compMS2", function(object, value) {
  
  object@spectralDB <- value
  return(object)
})

