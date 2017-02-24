#' customized PMID search function adapted from PubChemWordcloud package v 0.3.2
#' @param keys character vector of compound names to search pubmed with
#' @param n numeric maximum number of results to return
#' @export
PMIDsearch <- function(keys=NULL, n=1000){
  searchUrl <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="',
                    gsub(" ", "+", keys), paste0('"&retmax=', n), 
                    '&tool="compMS2Miner&email="edmandsw@berkeley.edu"')
  hlpURL <- RCurl::getURL(searchUrl, .opts=RCurl::curlOptions(followlocation=TRUE))
  doc <- XML::xmlTreeParse(hlpURL, asText = TRUE)
  IdlistHlp <- unlist(doc[["doc"]][["eSearchResult"]]["IdList"])
  if('QuotedPhraseNotFound' %in% unlist(doc[["doc"]][["eSearchResult"]][["WarningList"]])){
  Idlist <- 0
  } else {
  Count <- unlist(doc[["doc"]][["eSearchResult"]][["Count"]])[3]
  Idlist <- c(Count, IdlistHlp[grep("value$", names(IdlistHlp))])
  }
  return(Idlist)
} # end function
