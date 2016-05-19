#' Publish your compMS2explorer app on shinyapps.io
#' 
#' you must hold a shinyapps.io account and have a token.
#' use this command from the package rsconnect to set up your token
#' \code{setAccountInfo(name='your_shinyapps.io_user_name', token='your_token', secret='your_secret')}
#' Following this the app can be deployed on your account
#' @param object a compMS2 class object 
#' @param appName character name for your new app (e.g. 'compMS2example').
#' @param writeDir character full path to a directory to save the compMS2 results
#' and the current shiny app to. This zip file can then be shared with others.
#' @param ... further arguments to the \code{\link{deployApp}} function
#' @export
setGeneric("publishApp", function(object, appName=NULL, writeDir=NULL, ...) standardGeneric("publishApp"))

setMethod("publishApp", signature = "CompMS2", function(object, appName=NULL, ...){
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
    if (appDir == ""){
      stop('Could not find example directory. Try re-installing "CompMS2miner".', call. = FALSE)
    }

    save(object, file=paste0(appDir, '/compMS2object.RData'))  
    
    if(!is.null(writeDir)){
    # create directory to bundle app in
    appDirWrite <- paste0(writeDir, '/', appName)
    dir.create(appDirWrite)
    # id shiny files and .RData
    shinyFilesTmp <- list.files(appDir, full.names = T)
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
    # deploy app
    shinyapps::deployApp(appDir = appDir, appName=appName, ...)
    }
   }
}) # end function
    