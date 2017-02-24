#' lipid abbreviations table for pubmed text mining
#' 
#' This dataset contains fields from the LMSD database classifications and 
#' abbreviations (\url{http://www.lipidmaps.org/data/classification/lipid_cns.html})
#' The variables are as follows:
#'
#' \enumerate{
#'  \item Class lipid class name this string is intended to be searched in PubMed.
#'  \item Abbreviation abbreviation for lipid class.
#'  \item regexpr an R regular expression to try to detect the lipid class in a
#'  database compound name. e.g. Searching PubMed for the string "SM(18:1/14:0)"
#'  for example will return no PubMed ids however searching using "sphingomyelin"
#'  will return a more representative number of PubMed abstract Ids. 
#'  }
#' @docType data
#' @keywords datasets
#' @name lipidAbbrev
#' @usage data(lipidAbbrev)
#' @source \url{http://www.lipidmaps.org/data/classification/lipid_cns.html} 
#' @format A data frame with 18 rows and 3 columns
NULL
# common phospholipid abbreviations
# lipidmaps
# dated 07/15/2016
# require(XML)
# htmlTables <- XML::readHTMLTable('http://www.lipidmaps.org/data/classification/lipid_cns.html', stringsAsFactors = FALSE)
# lipidAbbrev <- data.frame(stringsAsFactors = FALSE)
# for(i in 2:4){
# lipidTableTmp <- htmlTables[[i]]
# # Glycosphingolipids
# lipidTableTmp <- lipidTableTmp[grepl('Glycosphingolipids', lipidTableTmp$Class) == FALSE, , drop=FALSE]
# lysoIndx <- grepl('lyso', lipidTableTmp$Abbreviation)
# indxTmp <-  grepl('[A-Z]', lipidTableTmp$Examples) | lysoIndx
# lipidTableTmp <- lipidTableTmp[indxTmp, , drop=FALSE]
# lysoIndx <- grepl('lyso', lipidTableTmp$Abbreviation)
# 
# lipidTableTmp$Class <- gsub('glycerophospho', 'phosphatidyl', lipidTableTmp$Class, ignore.case = TRUE)
# lipidTableTmp$Class <- gsub('radylglycerolipids', 'glycerides', lipidTableTmp$Class)
# lipidTableTmp$Abbreviation <- gsub(' \\(.+', '', lipidTableTmp$Abbreviation)
# lipidTableTmp$Examples <- NULL
# lipidTableTmp$regexpr <- gsub('$', '\\\\(', lipidTableTmp$Abbreviation)
# lipidTableTmp$regexpr <- gsub('^', '\\^', lipidTableTmp$regexpr)
# if(any(lysoIndx)){
# lysoLipids <- lipidTableTmp[lysoIndx, , drop=FALSE]
# lysoLipids$Class <- paste0('lyso', lysoLipids$Class)
# lysoLipids$Abbreviation <- paste(paste0('L', lysoLipids$Abbreviation), 
#                                  paste0('lyso', lysoLipids$Abbreviation), sep=' ')
# lysoLipids$regexpr <- gsub(' ', '\\\\(|\\^',  lysoLipids$Abbreviation)
# lysoLipids$regexpr <- gsub('$', '\\\\(', lysoLipids$regexpr)
# lysoLipids$regexpr <- gsub('^', '\\^', lysoLipids$regexpr)
# 
# lipidTableTmp <- rbind(lipidTableTmp, lysoLipids)  
# }
# clIndx <- grep('^CL$', lipidTableTmp$Abbreviation)
# if(length(clIndx) > 0){
# lipidTableTmp$Class[clIndx] <- 'cardiolipin'  
# }
# lipidAbbrev <- rbind(lipidAbbrev, lipidTableTmp)
# }
# lipidAbbrev$Class <- tolower(gsub('s$', '', lipidAbbrev$Class))
