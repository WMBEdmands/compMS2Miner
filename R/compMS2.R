#' Constructor for compMS2 class object from a peak table and MS2 mzXML file(s)
#'
#' @description Matches MS1 features to MS2 spectra (.mzXML) files based on a 
#' mass-to-charge and retention time tolerance. Composite spectra and other data
#' can subsequently be visualized during any stage of the CompMS2miner
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
#'  If argument is not supplied a GUI (tcltk) file selection window will open and a .csv file can then be selected. 
#' @param nSlaves numeric Number of cores for parallel computation.
#' @param mode character Ionisation polarity must be either 'pos' or 'neg'.
#' @param precursorPpm numeric Parts per million mass accuracy to match MS1 features to MS2 spectra (ppm) 
#' @param ret numeric retention time tolerance to match MS1 features to MS2 spectra (+/- seconds). 
#' @param TICfilter numeric Minimum Total Ion Current to consider an MS2 spectrum. Any MS2 scan
#' below this threshold will not be considered. 
#' 
#' @return A compMS2 object                   
#' @export
compMS2 <-  function(MS1features = NULL, mzXMLdir = NULL,
                     nSlaves = NULL, mode = "pos", 
                     precursorPpm = 10, ret = 10, TICfilter = 10000){
  
  message("creating compMS2 object in ", ifelse(mode == "pos","positive", 
                                                "negative"), 
          " ionisation mode")
  flush.console()
  # if foreach package not installed
  if(!require(foreach)){
    install.packages('foreach')
    if(!require(foreach)){
      stop('Unable to install the foreach package which is required for correct functioning of the compMS2miner package...\n')
    }
  }
  # if foreach package not installed
  if(!require(Rcpp)){
    install.packages('Rcpp')
    if(!require(Rcpp)){
      stop('Unable to install the Rcpp package which is required for correct functioning of the compMS2miner package...\n')
    }
    }
  
  # install shiny if necessary
  if(!require(shiny)){
    install.packages('shiny')
    if(!require(shiny)){
      stop('Unable to install the rCharts package which is required for correct functioning of the compMS2explorer application...\n')  
    }
  }
  # install rCharts/ devtools if necessary
  # if(!require(rCharts)){
  #   if(!require(devtools)){
  #     install.packages('devtools')
  #   }
  #   if(!require(devtools)){
  #     stop('Unable to install the devtools package which is required for correct functioning of the compMS2explorer application...\n')    
  #   }
  #   devtools::install_github('ramnathv/rCharts')
  #   if(!require(rCharts)){
  #     stop('Unable to install the rCharts package which is required for correct functioning of the compMS2explorer application...\n')  
  #   }  
  # } 
  # set proxy settings
#   setInternet2(TRUE) 
  # set global options
  options(stringsAsFactors = F)
  # if is.null mzXMLdir select mzXML file containing directory
  if(is.null(mzXMLdir)){
  message("Select your .mzXML data directory")
  flush.console()
  mzXMLdir <- tcltk::tk_choose.dir(default = "",  
                                    caption = "1. Select your .mzXML data directory")
  }
  # if is.null MS1features select MS1 feature table file
  if(is.null(MS1features)){
    message("Select your MS1 feature (.csv) file")
    flush.console()
    
  MS1features <- tcltk::tclvalue(tcltk::tkgetOpenFile(
       initialdir = mzXMLdir,
        filetypes = "{{Comma delimited} {.csv}} {{All files} *}", 
            title = "2. select the MS1 features (.csv) file you wish to match"))
  }
  
  if(is.character(MS1features)){
  # read in MS1features
  message("Reading MS1 feature table...")
  flush.console()
  MS1features <- as.data.frame(data.table::fread(MS1features, sep=",", 
                                                   header=T, stringsAsFactors=F))
  message("...Done")
  flush.console()
  }
  
  # sort dataframe by unique ID/ EIC number
  MS1features <- MS1features[order(MS1features[, 1]), ]
  row.names(MS1features) <- seq(1, nrow(MS1features), 1)
  
  # identify all mzXML files in raw-data directory
  MS2files <- list.files(path=mzXMLdir, pattern = "*.mzXML$", full.names=T)
  message(paste0(length(MS2files), 
         " MS2 (.mzXML) files were detected within the directory..."))
  flush.console()

  if(!is.null(nSlaves)){
    
    message(paste0("Starting SNOW cluster with ", nSlaves, " local sockets..."))
    flush.console()
    cl <- parallel::makeCluster(nSlaves) 
    doSNOW::registerDoSNOW(cl)
    
    message("matching MS1 peak table features to the following MS2 files: ")
    flush.console()
    mS2message <- sapply(basename(MS2files), function(x){ 
      message(x)
      flush.console()})
    # foreach and dopar from foreach package
    Results <- foreach(i=1:length(MS2files), .packages=c('mzR')) %dopar% {
      compMS2create(MS2file=MS2files[i], MS1features=MS1features,
                    TICfilter=TICfilter,  precursorPpm=precursorPpm, 
                    ret=ret)
    }
    # stop SNOW cluster
    parallel::stopCluster(cl) 
  } else {
    # create list to store results
    Results <- vector("list", length(MS2files))
    
    for(i in 1:length(MS2files)){
    Res.tmp <- compMS2create(MS2file=MS2files[i], TICfilter=TICfilter, 
                             MS1features=MS1features,
                             precursorPpm=precursorPpm, ret=ret)
    Results[[i]] <- Res.tmp
    }
  }
  
  Results <- unlist(Results, recursive = F)
  compMS2 <- new("CompMS2")
  filePaths(compMS2) <- MS2files
  compSpectra(compMS2) <- lapply(Results, function(x) x$spectra)
  metaData(compMS2) <- lapply(Results, function(x) x$metaData)
  # MS1 feature table
  # MS1features(compMS2) <- MS1features
  Parameters(compMS2) <- data.frame(nSlaves=ifelse(is.null(nSlaves), 0, nSlaves),
                                    mode=mode, precursorPpm=precursorPpm,
                                    ret=ret, TICfilter=TICfilter)
  return(compMS2)
  } # end CompMS2obj function

setMethod("show", "CompMS2", function(object) {
  if(length(object@file.paths) > 0){
    cat("A \"CompMS2\" class object derived from", length(object@file.paths), 
        "MS2 files \n\n")

    Acc.tmp <- sapply(object@metaData, function(x){ 
      c(mean(abs(unlist(x[grep("ppmDiff", names(x))]))),
        mean(abs(unlist(x[grep("rtDiff", names(x))]))),
      length(unlist(x[grep("rtDiff", names(x))])))})
    ionsCount <- sum(sapply(object@compSpectra, function(x) nrow(x)))
    
    spec.names <- names(object@compSpectra)
    cat(length(unique(gsub(".+_", "", spec.names))), 
        "MS1 features were matched to", sum(Acc.tmp[3, ]), "MS2 precursor scans\n")
    cat("containing",ionsCount, "ion features\n\n")    
    
    cat("Average ppm match accuracy:",
        round(mean(Acc.tmp[1, ]), digits = 3) ,"\n")
    cat("with a ppm mass accuracy tolerance of (+/-)", object@Parameters$precursorPpm,
        "\n\n")
    
    cat("Average retention time match accuracy:",
        round(mean(Acc.tmp[2, ]), digits = 2) ,"seconds\n")
    cat("with a retention time tolerance of (+/-)", object@Parameters$ret,
        " seconds\n\n")
    
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  } else {
    cat("A new empty\"CompMS2\" class object")
  }
})


# split method e.g. cut(1:length(x@compSpectra), 2)
setMethod("split", "CompMS2", function(x, f){
stopifnot(class(x) == 'CompMS2')
  
  splitX <- vector('list', nlevels(f))
  for(i in 1:nlevels(f)){
  tmpCompMS2 <- new('CompMS2')
  indxTmp <- which(f == levels(f)[i])
  tmpCompMS2@compSpectra <- x@compSpectra[indxTmp]
  tmpCompMS2@metaData <- x@metaData
  # tmpCompMS2@MS1features <- x@MS1features
  tmpCompMS2@DBanno <- x@DBanno[indxTmp]
  tmpCompMS2@BestAnno <- x@BestAnno[indxTmp]
  subStrIndx <- x@subStrAnno$compSpecName %in% names(tmpCompMS2@compSpectra)[indxTmp]
  tmpCompMS2@subStrAnno <- x@subStrAnno[subStrIndx, , drop=F]
  tmpCompMS2@Comments<- x@Comments[indxTmp]
  tmpCompMS2@file.paths <- x@file.paths
  tmpCompMS2@Parameters <- x@Parameters
  tmpCompMS2@couchDBconn <- x@couchDBconn
  tmpCompMS2@MetFrag <- x@MetFrag[indxTmp]
  tmpCompMS2@CFModelling <- x@CFModelling[indxTmp]
  splitX[[i]] <- tmpCompMS2
  names(splitX)[i] <- i
  }
  return(splitX)
})

setGeneric("filePaths", function(object, ...) standardGeneric("filePaths"))

setMethod("filePaths", "CompMS2", function(object) object@file.paths)

setGeneric("filePaths<-", function(object, value) standardGeneric("filePaths<-"))

setReplaceMethod("filePaths", "CompMS2", function(object, value) {
  object@file.paths <- value
  return(object)
})


setGeneric("compSpectra", function(object, ...) standardGeneric("compSpectra"))

setMethod("compSpectra", "CompMS2", function(object) object@compSpectra)

setGeneric("compSpectra<-", function(object, value) standardGeneric("compSpectra<-"))

setReplaceMethod("compSpectra", "CompMS2", function(object, value) {
  object@compSpectra <- value
  return(object)
})

setGeneric("metaData", function(object, ...) standardGeneric("metaData"))

setMethod("metaData", "CompMS2", function(object) object@metaData)

setGeneric("metaData<-", function(object, value) standardGeneric("metaData<-"))

setReplaceMethod("metaData", "CompMS2", function(object, value) {
  object@metaData  <- value
  return(object)
})

# setGeneric("MS1features", function(object, ...) standardGeneric("MS1features"))
# 
# setMethod("MS1features", "CompMS2", function(object) object@MS1features)
# 
# setGeneric("MS1features<-", function(object, value) standardGeneric("MS1features<-"))
# 
# setReplaceMethod("MS1features", "CompMS2", function(object, value) {
#   object@MS1features <- value
#   return(object)
# })

setGeneric("DBanno", function(object, ...) standardGeneric("DBanno"))

setMethod("DBanno", "CompMS2", function(object) object@DBanno)

setGeneric("DBanno<-", function(object, value) standardGeneric("DBanno<-"))

setReplaceMethod("DBanno", "CompMS2", function(object, value) {
  object@DBanno <- value
  return(object)
})

setGeneric("BestAnno", function(object, ...) standardGeneric("BestAnno"))

setMethod("BestAnno", "CompMS2", function(object) object@BestAnno)

setGeneric("BestAnno<-", function(object, value) standardGeneric("BestAnno<-"))

setReplaceMethod("BestAnno", "CompMS2", function(object, value) {
  object@BestAnno <- value
  return(object)
})

setGeneric("subStrAnno", function(object, ...) standardGeneric("subStrAnno"))

setMethod("subStrAnno", "CompMS2", function(object) object@subStrAnno)

setGeneric("subStrAnno<-", function(object, value) standardGeneric("subStrAnno<-"))

setReplaceMethod("subStrAnno", "CompMS2", function(object, value) {
  object@subStrAnno <- value
  return(object)
})

setGeneric("MetFrag", function(object, ...) standardGeneric("MetFrag"))

setMethod("MetFrag", "CompMS2", function(object) object@MetFrag)

setGeneric("MetFrag<-", function(object, value) standardGeneric("MetFrag<-"))

setReplaceMethod("MetFrag", "CompMS2", function(object, value) {
  object@MetFrag <- value
  return(object)
})

setGeneric("CFModelling", function(object, ...) standardGeneric("CFModelling"))

setMethod("CFModelling", "CompMS2", function(object) object@CFModelling)

setGeneric("CFModelling<-", function(object, value) standardGeneric("CFModelling<-"))

setReplaceMethod("CFModelling", "CompMS2", function(object, value) {
  object@CFModelling <- value
  return(object)
})

setGeneric("Comments", function(object, ...) standardGeneric("Comments"))

setMethod("Comments", "CompMS2", function(object) object@Comments)

setGeneric("Comments<-", function(object, value) standardGeneric("Comments<-"))


setReplaceMethod("Comments", "CompMS2", function(object, value) {
  object@Comments <- value
  return(object)
})

setGeneric("Parameters", function(object, ...) standardGeneric("Parameters"))

setMethod("Parameters", "CompMS2", function(object) object@Parameters)

setGeneric("Parameters<-", function(object, value) standardGeneric("Parameters<-"))

setReplaceMethod("Parameters", "CompMS2", function(object, value) {
  
  object@Parameters <- value
  return(object)
})

setGeneric("couchDBconn", function(object, ...) standardGeneric("couchDBconn"))

setMethod("couchDBconn", "CompMS2", function(object) object@couchDBconn)

setGeneric("couchDBconn<-", function(object, value) standardGeneric("couchDBconn<-"))

setReplaceMethod("couchDBconn", "CompMS2", function(object, value) {
  
  object@couchDBconn <- value
  return(object)
})

