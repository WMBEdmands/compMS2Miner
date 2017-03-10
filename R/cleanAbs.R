#' Adapted from PubMedWordCloud (cleanAbstracts) to work with compMS2Miner
#' @param Abs	output of getAbs, or just a paragraph of text
#' @export
cleanAbs <- function(Abs, rmNum = TRUE, tolw = TRUE, toup = FALSE, 
          rmWords = TRUE, yrWords = NULL, stemDoc = FALSE){
  abstTxt <- tm::Corpus(tm::VectorSource(Abs))
  text2.corpus = tm::tm_map(abstTxt, tm::removePunctuation)
  if (rmNum == TRUE) {
    text2.corpus = tm::tm_map(text2.corpus, function(x) tm::removeNumbers(x))
  }
  if (tolw == TRUE) {
    text2.corpus = tm::tm_map(text2.corpus, tolower)
  }
  if (toup == TRUE) {
    text2.corpus = tm::tm_map(text2.corpus, toupper)
  }
  if (rmWords == TRUE) {
    text2.corpus = tm::tm_map(text2.corpus, tm::removeWords, tm::stopwords("english"))
    if (!is.null(yrWords)) {
      text2.corpus = tm::tm_map(text2.corpus, tm::removeWords, 
                            yrWords)
    }
  }
  if (stemDoc == TRUE) {
    text2.corpus = tm_map(text2.corpus, tm::stemDocument)
  }
  text2.corpus <- tm_map(text2.corpus, tm::PlainTextDocument)
  
  # tdm <- TermDocumentMatrix(text2.corpus)
  indWords <- gsub(' ', '', unlist(strsplit(text2.corpus$content$content, ' ')))
  indWords <- indWords[indWords != '']
  tdm <- table(indWords)
  v <- sort(rowSums(m), decreasing = TRUE)
  d <- data.frame(word = names(v), freq = v)
  return(d)
} # end function
