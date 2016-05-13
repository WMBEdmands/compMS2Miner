#' Publish your compMS2explorer app on shinyapps.io
#' 
#' you must hold a shinyapps.io account and have a token.
#' use this command from the package rsconnect to set up your token
#' \code{rsconnect::setAccountInfo(name='your_shinyapps.io_user_name', token='your_token', secret='your_secret')}
#' Following this the app can be deployed on your account
#' @param object a compMS2 class object 
#' @param appName a name for your new app.
#' @param ... further arguments to the \code{\link{shinyapps::deployApp}} function
#' @export
setGeneric("publishApp", function(object, appName=NULL, ...) standardGeneric("publishApp"))

setMethod("publishApp", signature = "CompMS2", function(object, appName=NULL, ...){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(object@file.paths) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    
    if(!require(rsconnect)){
      stop('The package rsconnect must be installed to send you app to shinyapps.io...')
    }
    if(!require(shinyapps)){
      stop('The package shinyapps must be installed to send you app to shinyapps.io...')
    }
    
    appDir <- system.file("shiny-apps", "compMS2explorer", package = "CompMS2miner")
    if (appDir == "") {
      stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
    }
    save(object, file=paste0(appDir, '/compMS2object.RData'))  
    # deploy app
    shinyapps::deployApp(appDir = appDir, appName=appName, ...)
   }
}) # end function
    