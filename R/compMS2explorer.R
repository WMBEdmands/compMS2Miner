#' Visualize your compMS2Miner results output using a shiny app.
#' 
#' @param object a compMS2 class object or a character full path to a compMS2Miner zip file
#' @param browserLaunch logical launch app in web browser (default = TRUE).
#' @export
setGeneric("compMS2Explorer", function(object, ...) standardGeneric("compMS2Explorer"))

setMethod("compMS2Explorer", signature = "character", function(object, 
                                                                browserLaunch = TRUE){
  if(file.exists(object)){
    outdir <- tempfile(pattern = "compMS2Miner")
    dir.create(outdir)
    
    pathTmp <- utils::unzip(object, exdir = outdir)
    shiny::runApp(dirname(pathTmp[1]), launch.browser = browserLaunch)
  } else {
    stop('a full path to a compMS2Miner .zip file archive must be supplied')
  }
})
    
setMethod("compMS2Explorer", signature = "compMS2", function(object, 
                                                                browserLaunch = TRUE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(filePaths(object)) == 0) {
    stop("The CompMS2 class file is empty")
  } else {
    if(!require(shiny)){
      stop('The package shiny must be installed to use the compMS2Explorer function...')
    }
    appDir <- system.file("shiny-apps", "compMS2Explorer", package = "compMS2Miner")
      if (appDir == "") {
        stop("Could not find example directory. Try re-installing `compMS2Miner`.", call. = FALSE)
      }
    # add readOnly = F to parameters
    Parameters(object)$readOnly <- FALSE
    # create temporary directory to create zip
    outDir <- tempfile(pattern = "compMS2Miner")
    dir.create(outDir)
    wwwOutDir <- paste0(outDir, '/www')
    dir.create(wwwOutDir)
    # copy latest version of shiny app from package
    filesMoved <- file.copy(dir(appDir, full.names = TRUE, pattern = '\\.R$'), outDir, overwrite = TRUE)
    filesMoved <- c(filesMoved, file.copy(dir(paste0(appDir, '/www'), full.names = TRUE, pattern = '\\.mp4$|\\.png$'), wwwOutDir, overwrite = TRUE))
    if(any(filesMoved == FALSE)){
      stop('The shiny-app file(s) were not copied to the bundle please check the compMS2Miner package is properly installed.\n')
    }
    
    save(object, file=paste0(outDir, '/compMS2object.RData'))  
      shiny::runApp(outDir, display.mode = "normal", launch.browser = browserLaunch)
  }
}) # end function
