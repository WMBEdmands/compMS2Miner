#' run github shiny modified from shiny and devtools
#' @param repo character github username and repository name. in the form "username/repositoryName"
#' @param subdir character sub-directory of the repo containing the shiny and data.
#' @param dirPath character full-path to a directory in which to save the contents of the zip file. If unsupplied shiny app will be opened from a temporary directory. 
#' @param auth_token character private repo authorization token. 
#' @param browserLaunch logical launch app in web browser (default = TRUE).
#'
#' @export
runGitHubApp <- function(repo=NULL, subdir=NULL, dirPath=NULL, auth_token=NULL, browserLaunch=TRUE){
  # error handling
  stopifnot(is.character(repo))
  
    res <- strsplit(repo, "/")[[1]]
    if(length(res) != 2){ 
      stop("'repo' must be of the form 'username/repo'")
    }
    username <- res[1]
    repo <- res[2]
  
  remote <- devtools:::remote("github", host = "api.github.com", repo = repo, subdir =
                              subdir, username = username, ref = NULL, sha = NULL, 
                              auth_token = auth_token)
  
  bundle <- devtools:::remote_download(remote, quiet = FALSE)
  on.exit(unlink(bundle), add = TRUE)
  outdir <- tempfile(pattern = "compMS2Miner")
  dir.create(outdir)
  
  pathTmp <- utils::unzip(bundle, exdir = outdir) 
  # if sub directory not supplied then print list of options to console
  if(is.null(subdir)){
   availOpts <- basename(pathTmp) 
   availOpts <- availOpts[grep('\\.zip$', availOpts)]
   availOpts <- gsub('\\.zip', '', availOpts)
   message('\n"subdir" argument not supplied. Available directories within the repo include:\n',
           paste0(availOpts, '\n'), '\nPlease type a directory name without quotations and press [enter] to continue:')
   flush.console()
   subdir <- readline()
  }
  tmpIndx <- grepl(paste0(subdir, '.zip$'), pathTmp)
  if(!any(tmpIndx)){
    stop(subdir, ' sub-directory name not found please check and try again...')
  }
  tmpAppDir <- tempfile(pattern="compMS2Miner")
  if(!is.null(dirPath)){
  dirPath <- gsub("/$|\\\\$", '', dirPath)
  }
  appPathTmp <- utils::unzip(pathTmp[tmpIndx], exdir=ifelse(is.null(dirPath), tmpAppDir, dirPath))
  object <- shiny::runApp(dirname(appPathTmp[1]), launch.browser = browserLaunch)
  return(object)
} # end function
