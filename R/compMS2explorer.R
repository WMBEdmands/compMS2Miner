#' Visualize your CompMS2miner results output using a shiny app.
#' 
#' @param object a compMS2 class object 
#' @param browserLaunch logical launch app in web browser (default = TRUE).
#' @export
setGeneric("compMS2explorer", function(object, ...) standardGeneric("compMS2explorer"))

setMethod("compMS2explorer", signature = "CompMS2", function(object, browserLaunch = TRUE){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(object@file.paths) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    if(!require(shiny)){
      stop('The package shiny must be installed to use the compMS2explorer function...')
    }
    appDir <- system.file("shiny-apps", "compMS2explorer", package = "CompMS2miner")
      if (appDir == "") {
        stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
      }
    # add readOnly = F to parameters
    Parameters(object)$readOnly <- FALSE
    # create temporary directory to create zip
    outDir <- tempfile(pattern = "CompMS2miner")
    dir.create(outDir)
    
    # copy latest version of shiny app from package
    filesMoved <- file.copy(paste0(appDir, c('/server.R', '/ui.R', '/global.R')), outDir, overwrite = T)
    if(any(filesMoved == F)){
      stop('The shiny-app file(s) were not copied to the bundle please check the CompMS2miner package is properly installed.\n')
    }
    
    save(object, file=paste0(outDir, '/compMS2object.RData'))  
      shiny::runApp(outDir, display.mode = "normal", launch.browser = browserLaunch)
    }
}) # end function