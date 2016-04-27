#' return cleaned abstracts from pubmed from searched key words
#' @param keys character vector of compound names to search pubmed with
#' @param n numeric maximum number of results to return
#' @param maxChar numeric maximum number of characters in cleaned abstract words to return.
#' @param ... further arguments to the \code{\link{cleanAbstracts}} function.
#' @return a list containing 3 named elements:
#' 1. titles character vector of Abstract title(s)
#' 2. abs character vector of abstract text(s).
#' 3. clAbs clean abstract word frequency data.frame with column names 'word' and 'freq'. 
#' @seealso \code{\link{PubMedWordcloud}}, \code{\link{getAbstracts}}, \code{\link{cleanAbstracts}}.
#' @export
pubMedSearch <- function(keys=NULL, n=1000, maxChar=50, ...){
  #error handling
  if(is.null(keys)){
    stop('argument keys is missing with no default')
  }
  message('searching pubMed for the following key words :\n', paste0(keys, "\n"))
  flush.console()
  # search key words against pub med abstracts 
  PMIDs <- PMIDsearch(keys, n)
  message(length(PMIDs), ' pubmed IDs returned')
  flush.console()
  
  if(length(PMIDs) > 0)
  { 
    message('obtaining abstract text and titles from pubmed...')
    flush.console()
    # obtain abstract text
    Abs <- PubMedWordcloud::getAbstracts(PMIDs)
    # return titles
    titles <- getTitles(PMIDs)
    if(length(Abs) > 0)
    {
      message('cleaning abstracts (removing punctuations, numbers, translate characters to lower or upper case, remove stopwords, stemming words...')
      flush.console()
      
      ClAbs <- PubMedWordcloud::cleanAbstracts(Abs, ...)
    
      # only keep word which are less than max characters
      ClAbs <- ClAbs[which(sapply(as.character(ClAbs$word), nchar) < maxChar), ]
      # return results
      return(list(titles=titles, Abstracts=Abs, ClAbs=ClAbs))
   
    } else {
      stop("No abstract text was returned for the keyword(s) :\n", paste0(keys, "\n")) 
    }
  } else {
  stop("No pubmed ids returned for the keyword(s) :\n", paste0(keys, "\n")) 
  }
} # end function