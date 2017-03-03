#' Publish your compMS2Explorer app on shinyapps.io
#' 
#' you must hold a shinyapps.io account and have a token.
#' use this command from the package rsconnect to set up your token
#' \code{setAccountInfo(name='your_shinyapps.io_user_name', token='your_token', secret='your_secret')}
#' Following this the app can be deployed on your account
#' @param object a compMS2 class object 
#' @param appName character name for your new app (e.g. 'compMS2Example').
#' @param writeDir character full path to a directory to save the compMS2 results
#' and the current shiny app zip file to. This zip file can then be shared with others.if this argument is not supplied the results will be deployed on shinyapps.io.
#' @param addFiles character vector of full paths to files which will be included
#' in the zip file or bundle to shinyapps.io. For example code used to generate 
#' compMS2 results. The default is to include at minimum a text file containing 
#' the output of \link{sessionInfo()} this is intended to maintain reproducibility
#' of published results. 
#' @param ... further arguments to the \code{\link{deployApp}} function
#' @export
setGeneric("publishApp", function(object, appName=NULL, writeDir=NULL, addFiles=NULL, ...) standardGeneric("publishApp"))

setMethod("publishApp", signature = "compMS2", function(object, appName=NULL,  writeDir=NULL, addFiles=NULL, ...){
  # error handling
  if(!require(splashR)){
  stop('Please install splashR from GitHub using devtools: devtools::install_github("berlinguyinca/spectra-hash", subdir="splashR") this will be used to generated a unique splash code (spectral hash codes) for each spectrum.')
}
  if(class(object) != "compMS2"){
    stop("argument object is not a CompMS2 class object")
  } else if (length(filePaths(object)) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    
    if(!require(rsconnect)){
      stop('The package rsconnect must be installed to send your app to shinyapps.io...')
    }
    # if(!require(shinyapps)){
    #   stop('The package shinyapps must be installed to send your app to shinyapps.io...')
    # }
    
    message('Generating splash codes and adding to compMS2 object.\n')
    flush.console()
    metaData(object) <- sapply(1:length(compSpectra(object)), function(x){
      splashCode <- splashR::getSplash(compSpectra(object)[[x]][, c('mass', 'intensity')])
      # add to metaData
      metaDataTmp <- metaData(object)[[x]]
      # existing splash
      exSplash <- grep('splash', metaDataTmp)
      if(length(exSplash) == 0){
      entName <- paste0(gsub('_([^_]*)$', '', names(metaDataTmp)[1]), '_splashCode')
      metaDataTmp <- c(metaDataTmp, splashCode)
      names(metaDataTmp)[length(metaDataTmp)] <- entName 
      } else {
        metaDataTmp[[exSplash]] <- splashCode
      }
      return(metaDataTmp)
      })
    names(metaData(object)) <- names(compSpectra(object))
    appDir <- system.file("shiny-apps", "compMS2Explorer", package = "compMS2Miner")
    if(appDir == ""){
      stop('Could not find shiny-apps directory. Try re-installing "compMS2Miner".', call. = FALSE)
    }
    
    # create temporary directory to create zip
    outDir <- tempfile(pattern = "compMS2Miner")
    dir.create(outDir)
    wwwOutDir <- paste0(outDir, '/www')
    dir.create(wwwOutDir)
    
    if(is.character(addFiles)){
    message('moving the following additional file(s) to the application bundle "www" directory:\n',
            paste0(basename(addFiles), '\n'))
    flush.console()
    filesMoved <- file.copy(addFiles, wwwOutDir, overwrite = TRUE)
    if(any(filesMoved == FALSE)){
    stop('The file(s):\n', paste0(addFiles[filesMoved == FALSE], '\n'), 'were not copied to the bundle please check the path is correct.\n')
    }
    }  
    # generate sessionInfo text file
    seshInfoTxt <- paste0(outDir, '/sessionInfo_', gsub('-', '', Sys.Date()),
                          '.txt')
    writeLines(capture.output(sessionInfo()), seshInfoTxt)
    
    # copy latest version of shiny app from package
    filesMoved <- file.copy(dir(appDir, full.names = TRUE, pattern = '\\.R$'), outDir, overwrite = TRUE)
    filesMoved <- c(filesMoved, file.copy(dir(paste0(appDir, '/www'), full.names = TRUE, pattern = '\\.mp4$|\\.png$'), wwwOutDir, overwrite = TRUE))
    if(any(filesMoved == FALSE)){
      stop('The shiny-app file(s) were not copied to the bundle please check the compMS2Miner package is properly installed.\n')
    }
    
    if(!is.null(writeDir)){
    # add readOnly = F to parameters when zipped so collaborators can create their own edits
    Parameters(object)$readOnly <- FALSE
    save(object, file=paste0(outDir, '/compMS2object.RData'))  
    
    # create directory to bundle app and create zip in
    appDirWrite <- paste0(writeDir, '/', appName)
    dir.create(appDirWrite)
    # id shiny files and .RData
    shinyFilesTmp <- list.files(outDir, full.names = TRUE)
    # remove the dreaded desktop.ini
    shinyFilesTmp <- shinyFilesTmp[grepl('desktop\\.ini$', shinyFilesTmp) == FALSE]
    destZipFile <- paste0(appDirWrite, '/', appName, '.zip')
    fileAlreadyExists <- file.exists(destZipFile)
    if(fileAlreadyExists == TRUE){
      message('destination .zip file:\n\n', destZipFile, '\n\nalready exists...\nDo you wish to overwrite?\n')
      flush.console()
      ovWrite <- readline(prompt="Y/N: ")
      if(grepl('Y', ovWrite, ignore.case = TRUE)){
        tmpWd <- getwd()
        setwd(dirname(shinyFilesTmp))
        fileZipLog <- zip(destZipFile, basename(shinyFilesTmp))
        setwd(tmpWd)
      } else {
        message('stopping...\n')
        flush.console()
        opt <- options(show.error.messages=FALSE) 
        on.exit(options(opt)) 
        stop() 
      }
    } else {
    tmpWd <- getwd()
    setwd(dirname(shinyFilesTmp))
    fileZipLog <- zip(destZipFile, basename(shinyFilesTmp))
    setwd(tmpWd)
    }
    
    # # packrat shiny app
    # packrat::init(appDirWrite) 
    # 
    } else {
    # add readOnly = TRUE to parameters when published so viewers cannot edit
    Parameters(object)$readOnly <- TRUE
    save(object, file=paste0(outDir, '/compMS2object.RData'))  
      
    # deploy app
    shinyapps::deployApp(appDir = outDir, appName=appName, ...)
    }
   }
}) # end function
    
