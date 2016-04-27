#' Search pubmed compound for molecular formula using the pubchem power user gateway (PUG)
#' @param mf molecular formula character vector of length one e.g. 'C10H21N'. 
#' @return returns a character vector of pubmed compound cids matching the molecular
#' formula
#' @export
mfSearchPUG <- function(mf='C10H21N'){
  qUrl <- paste0('http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/', mf, '/XML')
  parsedhtml <- XML::htmlParse(qUrl)
  
  listKeyTmp <- parsedhtml['//listkey', fun=xmlValue][[1]]
  
  checkUrl <- paste0('http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/', 
                     listKeyTmp, '/cids/XML')
  parsedCheck <- XML::htmlParse(checkUrl)
  
  messageTmp <- tryCatch(parsedCheck['//message', fun=xmlValue][[1]], error=function(cond){
    message('query complete')
  })
  
  message(paste0(messageTmp, '...\n'))
  flush.console()
  # loop until query completed
  while(messageTmp == 'Your request is running'){
    Sys.sleep(3)
    parsedCheck <- tryCatch(XML::htmlParse(checkUrl), error=function(cond){
      return('no results returned')
    })
    
    if(is.character(parsedCheck)){
      message(paste0(parsedCheck, '...\n'))
      flush.console()
      break
    }
    messageTmp <- tryCatch(parsedCheck['//message', fun=xmlValue][[1]], error=function(cond){
      return('query complete')
    })
    # when query complete break loop
    if(messageTmp == 'query complete'){
      message(paste0(messageTmp, '...\n'))
      flush.console()
      break
    }
  }
  
  if(!is.character(parsedCheck)){
  cidsTmp <- parsedCheck['//cid', fun=xmlValue]
  
  # n cids returned
  message(length(cidsTmp), ' pubChem compound ids returned...\n')
  flush.console()
  return(as.numeric(unlist(cidsTmp)))
  } 
} # end function