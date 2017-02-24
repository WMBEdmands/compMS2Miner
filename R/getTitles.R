#' get PubMed title function adapted from PubChemWordcloud package v 0.3.2
#' @param pmid pubmed id number
#' @export
getTitles <- function(pmid){
  if(length(pmid) > 0){
    eDDownload <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
    hlp1 <- paste(eDDownload, paste(pmid, collapse = ",", 
                                    sep = ""), sep = "")
    hlp2 <- paste(hlp1, "&rettype=Abstract", sep = "")
    hlpURL <- RCurl::getURL(hlp2, .opts=RCurl::curlOptions(followlocation=TRUE))
    testDoc <- XML::xmlTreeParse(hlpURL, useInternalNodes = TRUE)
    topFetch <- XML::xmlRoot(testDoc)
    Title <- XML::xpathSApply(topFetch, "//ArticleTitle", XML::xmlValue)
  } else {
    Title<-c("Zero", "Articles", "Found")
  }
  return(Title)
}
