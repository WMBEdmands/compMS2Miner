#' Publish your compMS2explorer app on shinyapps.io
#' 
#' you must hold a shinyapps.io account and have a token.
#' use this command from the package rsconnect to set up your token
#' \code{setAccountInfo(name='your_shinyapps.io_user_name', token='your_token', secret='your_secret')}
#' Following this the app can be deployed on your account
#' @param object a compMS2 class object 
#' @param appName character name for your new app (e.g. 'compMS2example').
#' @param writeDir character full path to a directory to save the compMS2 results
#' and the current shiny app zip file to. This zip file can then be shared with others.if this argument is not supplied the results will be deployed on shinyapps.io.
#' @param addFiles character vector of full paths to files which will be included
#' in the zip file or bundle to shinyapps.io. For example code used to generate 
#' compMS2 results.
#' @param ... further arguments to the \code{\link{deployApp}} function
#' @export
setGeneric("publishApp", function(object, appName=NULL, writeDir=NULL, addFiles=NULL, ...) standardGeneric("publishApp"))

setMethod("publishApp", signature = "CompMS2", function(object, appName=NULL,  writeDir=NULL, addFiles=NULL, ...){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not a CompMS2 class object")
  } else if (length(object@file.paths) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    
    if(!require(rsconnect)){
      stop('The package rsconnect must be installed to send your app to shinyapps.io...')
    }
    if(!require(shinyapps)){
      stop('The package shinyapps must be installed to send your app to shinyapps.io...')
    }
    
    appDir <- system.file("shiny-apps", "compMS2explorer", package = "CompMS2miner")
    if(appDir == ""){
      stop('Could not find shiny-apps directory. Try re-installing "CompMS2miner".', call. = FALSE)
    }
    
    # create temporary directory to create zip
    outDir <- tempfile(pattern = "CompMS2miner")
    dir.create(outDir)
    
    if(is.character(addFiles)){
    message('moving the following additional files to the application bundle:\n',
            paste0(basename(addFiles), '\n'))
    flush.console()
    filesMoved <- file.copy(addFiles, outDir, overwrite = T)
    if(any(filesMoved == F)){
    stop('The file(s):\n', paste0(addFiles[filesMoved == F], '\n'), 'were not copied to the bundle please check the path is correct.\n')
    }
    }  
    # copy latest version of shiny app from package
    filesMoved <- file.copy(paste0(appDir, c('/server.R', '/ui.R', '/global.R')), outDir, overwrite = T)
    if(any(filesMoved == F)){
      stop('The shiny-app file(s) were not copied to the bundle please check the CompMS2miner package is properly installed.\n')
    }
    
    if(!is.null(writeDir)){
    # add readOnly = F to parameters when zipped so collaborators can create their own edits
    Parameters(object)$readOnly <- FALSE
    save(object, file=paste0(outDir, '/compMS2object.RData'))  
    
    # create directory to bundle app and create zip in
    appDirWrite <- paste0(writeDir, '/', appName)
    dir.create(appDirWrite)
    # id shiny files and .RData
    shinyFilesTmp <- list.files(outDir, full.names = T)
    # remove the dreaded desktop.ini
    shinyFilesTmp <- shinyFilesTmp[grepl('desktop\\.ini$', shinyFilesTmp) == F]
    destZipFile <- paste0(appDirWrite, '/', appName, '.zip')
    fileAlreadyExists <- file.exists(destZipFile)
    if(fileAlreadyExists == T){
      message('destination .zip file:\n\n', destZipFile, '\n\nalready exists...\nDo you wish to overwrite?\n')
      flush.console()
      ovWrite <- readline(prompt="Y/N: ")
      if(grepl('Y', ovWrite, ignore.case = T)){
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
    