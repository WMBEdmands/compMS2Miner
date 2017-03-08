#' return cleaned abstracts from pubmed from searched key words
#' @param keys character vector of compound names to search pubmed with
#' @param n numeric maximum number of results to return. The maximum and default is 500.
#' @param maxChar numeric maximum number of characters in cleaned abstract words to return.
#' @param ... further arguments to the \code{\link{cleanAbstracts}} function.
#' @return a list containing 3 named elements:
#' 1. titles character vector of Abstract title(s)
#' 2. abs character vector of abstract text(s).
#' 3. clAbs clean abstract word frequency data.frame with column names 'word' and 'freq'. 
#' @seealso PubMedWordcloud, \code{\link{getAbstracts}}, \code{\link{cleanAbstracts}}.
#' @export
pubMedSearch <- function(keys=NULL, n=500, maxChar=50, ...){
  #error handling
  if(is.null(keys)){
    stop('argument keys is missing with no default')
  }
  if(n > 500){
    stop('The maximum PMID key length is 500')
  }
  message('searching pubMed for the following key words :\n', paste0(keys, "\n"))
  flush.console()
  # search key words against pub med abstracts 
  PMIDs <- PMIDsearch(keys, n)
  message(length(PMIDs) - 1, ' pubmed IDs returned')
  flush.console()
  
  if(length(PMIDs) > 0)
  { 
    message('obtaining abstract text and titles from pubmed...')
    flush.console()
    # obtain abstract text
    Abs <- getAbs(PMIDs[-1])
    # return titles
    titles <- getTitles(PMIDs[-1])
    if(length(Abs) > 0){
      message('cleaning abstracts (removing punctuations, numbers, translate characters to lower or upper case, remove stopwords, stemming words...')
      flush.console()
      
      ClAbs <- PubMedWordcloud::cleanAbstracts(Abs, ...)
    
      # only keep word which are less than max characters
      ClAbs <- ClAbs[which(sapply(as.character(ClAbs$word), nchar) < maxChar), , drop=FALSE]
      # return results
      return(list(titles=titles, Abstracts=Abs, ClAbs=ClAbs))
   
    } else {
      stop("No abstract text was returned for the keyword(s) :\n", paste0(keys, "\n")) 
    }
  } else {
  stop("No pubmed ids returned for the keyword(s) :\n", paste0(keys, "\n")) 
  }
} # end function
