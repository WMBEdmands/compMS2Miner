#' customized PMID search function adapted from PubChemWordcloud package v 0.3.2
#' @param keys character vector of compound names to search pubmed with
#' @param n numeric maximum number of results to return
PMIDsearch <- function(keys=NULL, n=1000){
  searchUrl <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="',
                    gsub(" ", "+", keys), paste0('"&retmax=', n))
  hlpURL <- RCurl::getURL(searchUrl)
  doc <- XML::xmlTreeParse(hlpURL, asText = TRUE)
  IdlistHlp <- unlist(doc[["doc"]][["eSearchResult"]]["IdList"])
  Count <- unlist(doc[["doc"]][["eSearchResult"]][["Count"]])[3]
  Idlist <- c(Count, IdlistHlp[grep("value$", names(IdlistHlp))])
  return(Idlist)
} # end function