#' run github shiny modified from shiny and devtools
#' @param repo character github username and repository name. in the form "username/repositoryName"
#' @param subdir character sub-directory of the repo containing the shiny and data.
#' @param auth_token character private repo authorization token. 
#' @export
runGitHubApp <- function(repo=NULL, subdir=NULL, auth_token=NULL){
  # error handling
  stopifnot(is.character(repo))
  stopifnot(is.character(subdir))
  
    res <- strsplit(repo, "/")[[1]]
    if(length(res) != 2){ 
      stop("'repo' must be of the form 'username/repo'")
    }
    username <- res[1]
    repo <- res[2]
  
  remote <- devtools:::remote("github", host = "api.github.com", repo = repo, subdir =
                              subdir, username = username, ref = NULL, sha = NULL, 
                              auth_token = auth_token)
  
  bundle <- devtools:::remote_download(remote, quiet = F)
  on.exit(unlink(bundle), add = TRUE)
  outdir <- tempfile(pattern = "CompMS2miner")
  dir.create(outdir)
  
  pathTmp <- utils::unzip(bundle, exdir = outdir) 
  
  tmpIndx <- grepl(paste0('/', subdir, '/'), pathTmp)
  if(!any(tmpIndx)){
    stop(subdir, ' sub-directory name not found please check and try again...')
  }
  tmpAppDir <- tempfile(pattern='CompMS2miner')
  appPathTmp <- utils::unzip(pathTmp[tmpIndx], exdir=tmpAppDir)
  shiny::runApp(dirname(appPathTmp[1]))
} # end function