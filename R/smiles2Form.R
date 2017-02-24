#' Convert SMILES code to atomic formula
#' @param SMILES character vector of SMILES codes to convert
#' @return a character vector the formula(e).
smiles2Form <- function(SMILES=NULL){
  if(is.null(SMILES)){
    stop('SMILES argument is missing with no default')
  }
  if(!require(ChemmineR)){
    stop('ChemmineR package must be installed to use this function.')
  }
  data("exactMassEle")
  # identify replicates
  constEle <- gsub("[^[:alnum:] ]|[0-9]", "", SMILES)
  constEle <- gsub('([[:upper:]])', ' \\1', constEle)
  constEle <- gsub("([[:lower:]])([[:lower:]][[:lower:]])", "\\1 \\2", constEle)
  potDupForm <- sapply(strsplit(constEle, ' '), function(x){
    x <- x[x != '']
    tmpIdx <- x %in% exactMassEle$eleSymbol
    if(any(tmpIdx)){
    x <- c(unlist(strsplit(toupper(x[tmpIdx == FALSE]), '')), x[tmpIdx])
    }
   formNoImplH <- table(x)
   formNoImplH <- formNoImplH[order(names(formNoImplH))]
   formNoImplH <- paste(paste0(names(formNoImplH), formNoImplH), collapse = '')
   return(formNoImplH)
  })
  # in spite of errors only convert non duplicates
  nonRedSmiIdx <- duplicated(potDupForm) == FALSE
  sdfTmp <- suppressWarnings(ChemmineR::smiles2sdf(SMILES[nonRedSmiIdx]))
  formTmp <- ChemmineR::MF(sdfTmp, addH=TRUE)
  names(formTmp) <- potDupForm[nonRedSmiIdx]
  convFormulae <- formTmp[potDupForm]
  names(convFormulae) <- NULL
  return(convFormulae)
} # end function
