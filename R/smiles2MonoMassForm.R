#' Convert SMILES code to monoisotopic mass and formula
#' @param SMILES character vector of SMILES codes to convert
#' @return a named numeric vector of same length as the SMILES input containing the
#' monoisotopic mass(es) and named using the formula(e).
smiles2MonoMassForm <- function(SMILES=NULL){
  if(is.null(SMILES)){
    stop('SMILES argument is missing with no default')
  }
  if(!require(ChemmineR)){
    stop('ChemmineR package must be installed to use this function.')
  }
    sdfTmp <- suppressWarnings(ChemmineR::smiles2sdf(SMILES))
    formTmp <- ChemmineR::MF(sdfTmp, addH=TRUE)
    # monoisotopic mass
    eleGroupsStr <- toupper(formTmp)
    eleGroupsStr <- gsub('([[:upper:]])', ' \\1', formTmp)
    eleGroupsStr <- gsub('$', ' ', eleGroupsStr)
    eleGroupsStr <- gsub("^ ", '', eleGroupsStr)

    eleOnly <- gsub('[[:punct:]]|[0-9]', '', eleGroupsStr)
    nEleOnly <- gsub('[[:punct:]]|[A-z]', "", eleGroupsStr)
    eleOnly <- strsplit(eleOnly, ' ')
    nEleOnly <- lapply(strsplit(nEleOnly, ' '), function(x) as.numeric(ifelse(x == '', 1, x)))
    
    data("exactMassEle")
    massMostAbund <- apply(exactMassEle, 1, function(x){
      maxAbundTmp <- which.max(as.numeric(strsplit(x['natAbund'], ' ')[[1]]))
      mass <- as.numeric(strsplit(x['monoMass'], ' ')[[1]])[maxAbundTmp]
      return(mass)
    })
    names(massMostAbund) <- exactMassEle$eleSymbol
    # monoisotopic mass
    monoMass <- vector('numeric', length(eleOnly))
    for(i in 1:length(eleOnly)){
      eleOnly[[i]] <- unlist(mapply(rep, eleOnly[[i]], each=nEleOnly[[i]]))
      monoMass[i] <- sum(massMostAbund[eleOnly[[i]]])
    }
    names(monoMass) <- formTmp
    return(monoMass)
} # end function
