#' Adapted from PubMedWordCloud to work with compMS2Miner
#' @param PMID character vector of pubMed ids to get abstracts for.
#' @details if the query sequence is too long than 500 this function will not work
#' @export
getAbs <- function(PMID){
  if(!require(XML)){
    stop('The XML package is required to use this function')
  }
  if(length(PMID) > 500){
    stop('The maximum PMID length is 500')
  }
  if (length(PMID) > 0) {
    eDDownload <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
    hlp1 <- paste(eDDownload, paste(PMID, collapse = ",", 
                                    sep = ""), sep = "")
    hlp2 <- paste(hlp1, "&rettype=abstract", sep = "")
    hlpURL <- RCurl::getURL(hlp2, .opts=RCurl::curlOptions(followlocation=TRUE))
    testDoc <- XML::xmlTreeParse(hlpURL, useInternalNodes = TRUE)
    topFetch <- XML::xmlRoot(testDoc)
    abst <- XML::xpathSApply(topFetch, "//Abstract", xmlValue)
  }
  else {
    abst = c("Zero", "Articles", "Found")
  }
  return(abst)
} # end function
