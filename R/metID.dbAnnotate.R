#' Annotate unknown features in a compMS2 class object to database 
#' entries based on monoisotopic mass
#' 
#' @description unknown metabolite identification. MS1 features within a compMS2
#' object are matched against a metabolite database based on a mass accuracy
#' tolerance (ppm).
#' Possible electrospray adducts/ in-source fragments and substructure mass 
#' shifts are calculated for all data base entries. Substructure mass shifts 
#' are supplied as a named numeric vector or default is all the mass shifts in 
#' the internal substructure masses table (see default Substructure_masses).
#' Custom ESI adduct/ in-source fragments can be supplied as a character string
#' of names or the output of the function \link{adduct2mass}. See examples.
#' 
#' @param object a compMS2 class object.
#' @param ppm numeric mass accuracy error (ppm) to match MS1 features to data 
#' base entries. 
#' @param esiAdducts data.frame or character vector of custom electrospray 
#' adducts can be supplied as a character vector of electrospray adducts or as
#' a data.frame output from the function \link{adduct2mass}. e.g. 
#' c("[M-H-NH3-CO-COCH2-C4H6O]-", "[4M-H+Cl]2-", "[2M-H]-", "[3M-3H+Fe2+]-", 
#' "[M-H-CH2O]-", "[3M-2H]2-", "[M-H-CO2-C3H6]-", "[M-H+CH3COOH]-", "[3M-H]-", 
#' "[M-2H]2-")
#' @param metDB a metabolite data base (see default ?HMDB for table format).
#' Other currently available data bases include LMSD, DrugBank, T3DB and ReSpect 
#' databases. Matching using multiple databases is also possible. 
#' @param SubStrs named numeric vector of substructure mass shifts to consider. 
#' If argument is not supplied then no substructure mass shifts will be considered.
#' If the character "All" is supplied then all the substructure mass shifts 
#' contained in the internal table will be considered 
#' (see default \link{Substructure_masses} for more information). 
#' @param featureSubSet character vector of composite spectra names (e.g. CC_1, CC_2 etc.) otherwise the default is to perform database annotation on all composite spectra.  
#' @param includeElements character vector of element symbols (case-sensitive) 
#' to include in the metDB argument. Any structure containing an element which 
#' is not in this inclusion list will not be considered. see ?exactMassEle for 
#' the internal element table.
#' @param mixtures logical (default = FALSE), should mixtures be considered. Any
#' SMILES codes consisting of multiple non-covalently linked structures will not
#' be considered (e.g. salts)
#' @param MS1adducts logical (default = FALSE), should the adducts included in
#' the 4th column of the MS1feature table be used to reduce false positive 
#' assignments. The adducts identified will be used to guide the annotation process,
#' only the unique adducts/fragments will be used to calculate the expected
#' mass shift values. If a spectrum has no adduct identified then all of the
#' unique adducts identified in the dataset will be used. 
#' 
#' @return a compMS2 class object containing potential metabolite annotations.
#' @examples 
#' SubStrMassShift <- c(42.010565, 119.004101, 176.03209, 255.988909, 305.068159, 
#'                      57.021464, 161.014666, 79.956817)
#' names(SubStrMassShift) <- c("acetyl", "cysteine", "glucuronide", 
#'                             "glucuronide sulfate", "glutathione", "glycine",
#'                             "mercapturate", "sulfate")
#' # custom ESI adducts (default is to consider all ESI adducts/ in-source 
#' # fragments from supplementary material Beyond Profiling manuscript 
#' # see references for details)
#' # The function adduct2mass can interpret ESI adduct names and generate a 
#' # data.frame of expected mass shifts
#' customEsiAdducts <- c("[M-H-NH3-CO-COCH2-C4H6O]-", "[4M-H+Cl]2-", "[2M-H]-", 
#'                       "[3M-3H+Fe2+]-", "[M-H-CH2O]-", "[3M-2H]2-", 
#'                       "[M-H-CO2-C3H6]-", "[M-H+CH3COOH]-", "[3M-H]-", 
#'                       "[M-2H]2-")
#' 
#' compMS2Example <- metID(compMS2Example, "dbAnnotate", SubStrs=SubStrMassShift, 
#'                         esiAdducts=customEsiAdducts)
#'                                       
#' @references 
#' \enumerate{
#' \item Stanstrup, J., Gerlich, M., Dragsted, L.O. et al. 
#' Anal Bioanal Chem (2013) 405: 5037. doi:10.1007/s00216-013-6954-6
#' }
#' @seealso \link{adduct2mass}.
#' @export
setGeneric("metID.dbAnnotate", function(object, ...) standardGeneric("metID.dbAnnotate"))

setMethod("metID.dbAnnotate", signature = "compMS2", function(object, 
                                                              esiAdducts=NULL, 
                                                              ppm = NULL, 
                                                              metDB = HMDB, 
                                                              SubStrs = NULL,
                                                              featureSubSet=NULL,
                                                              includeElements=NULL,
                                                              mixtures=FALSE,
                                                              MS1adducts=FALSE){
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an compMS2 class object")
  }  
    if(is.null(ppm)){
      ppm <- Parameters(object)$precursorPpm
    }
  # add any MD or chemFP columns if necessary
  # availDbs <- c('HMDB', 'LMSD', 'DrugBank', 'T3DB', 'ReSpect')
  # includeMDChemFP <- TRUE
  # if(deparse(substitute(metDB))  %in% availDbs){
  #   includeMDChemFP <- FALSE
  # } 
  
    if(!all(c("Unique_DB_ID", "monoisotopic_weight", "name", "SMILES") %in% colnames(metDB))){
      stop('metDB colnames for unique database IDs, monoisotopic weights, names and canonical SMILES of data base
           entries must be named "Unique_DB_ID", "monoisotopic_weight", "name" and "SMILES" respectively' )
    }
  # remove non-covalent mixtures
  if(mixtures == FALSE){
    # remove non-covalently bound groups (e.g. mixtures and salts)
    nonCovIndx <- grepl('\\..+', metDB$SMILES)
    if(any(nonCovIndx)){
      dispIndx <- ifelse(sum(nonCovIndx) > 10, 10, sum(nonCovIndx))
      message('The following ', round(sum(nonCovIndx), 0), ' SMILES code(s) contain non-covalently bound groups and will not be considered:\n', 
              paste0(which(nonCovIndx)[1:dispIndx], '. ',  
                    (metDB$SMILES[nonCovIndx])[1:dispIndx], collapse='\n'), 
              ifelse(sum(nonCovIndx) > 10, '...(limited to first 10)\n', '\n'))
      flush.console()
      metDB <- metDB[nonCovIndx == FALSE, , drop=FALSE]
    }
  }
  
  # if necessary element inclusion list
  if(!is.null(includeElements)){
  nonStEle <- setdiff(exactMassEle$eleSymbol, includeElements)
  nonStEleRegEx <- gsub('^', ' ', paste0(nonStEle, collapse = ' | '))
  nonStEleRegEx <- gsub('$', ' ', nonStEleRegEx)
  # non-include elements
  # remove all punctuation character symbols and numbers
  eleGroupsStr <- gsub('[[:punct:]]|[0-9]', '', metDB$molecular_formula)
  # split into elements
  eleGroupsStr <- gsub('([[:upper:]])', ' \\1', eleGroupsStr)
  eleGroupsStr <- gsub('$', ' ', eleGroupsStr)
  nonInclIndx <- grepl(nonStEleRegEx, eleGroupsStr)
  if(any(nonInclIndx)){
    dispIndx <- ifelse(sum(nonInclIndx) > 10, 10, sum(nonInclIndx))
    message('\n\nThe following ', round(sum(nonInclIndx), 0), ' SMILES code(s) contain elements not present in the inclusion list and will not be considered:\n', 
            paste0(which(nonInclIndx)[1:dispIndx], '. ',  
                   (metDB$SMILES[nonInclIndx])[1:dispIndx], collapse='\n'), 
                   ifelse(sum(nonInclIndx) > 10, '...(limited to first 10)\n', 
                          '\n'))
    flush.console()  
    metDB <- metDB[nonInclIndx == FALSE, , drop=FALSE]
  }
  }
  # less than 3 character smiles remove
  elementsIndx <- nchar(metDB$SMILES) > 3
  if(any(elementsIndx == FALSE)){
    dispIndx <- ifelse(sum(elementsIndx == FALSE) > 10, 10, sum(elementsIndx == FALSE))
    message('\n\nThe following ', round(sum(elementsIndx == FALSE), 0), 
            ' SMILES code(s) are less than 4 characters in length and will not be considered:\n', 
            paste0(which(elementsIndx == FALSE)[1:dispIndx], '. ',  
                   (metDB$SMILES[elementsIndx == FALSE])[1:dispIndx], collapse='\n'), 
            ifelse(sum(elementsIndx == FALSE) > 10, '...(limited to first 10)\n', 
                   '\n'))
    flush.console()  
    metDB <- metDB[elementsIndx, , drop=FALSE]  
  }
  
  esiAdductsBySpec <- NULL  
  if(MS1adducts == TRUE){
  esiAdductsBySpec <- sapply(metaData(object), function(x){ 
    esiTmp <- unique(unlist(x[grep('_MS1_adduct$', names(x))]))
    esiTmp <- strsplit(esiTmp, ' ')[[1]]
    esiTmp <- esiTmp[grepl('\\[.*\\]', esiTmp)]
    if(!is.null(esiAdducts)){
      esiTmp <- c(esiTmp, esiAdducts)
    }
    esiTmp <- paste0(esiTmp, collapse = ' ')
    })
  names(esiAdductsBySpec) <- names(metaData(object))
  # unique adducts
  esiAdducts <- unique(unlist(strsplit(esiAdductsBySpec, ' ')))
  if(length(esiAdducts) == 0){
    stop('No MS1adducts were identified.')
  } else {
    message(length(esiAdducts), ' unique adducts/fragments were identified in the MS1-level data.\n')
    flush.console()
  }
  }
    if(!is.null(esiAdducts)){
      if(is.character(esiAdducts)){
        indxTmp <- grepl('\\[.*\\]', esiAdducts)
        if(!all(indxTmp)){
          stop('The following esiAdducts are not in the correct format:\n',
               paste0(esiAdducts[indxTmp == F], '\n'), '\ne.g. [4M-H+Cl]2-\n')
        }
        esiAdducts <- adduct2mass(esiAdducts)
      }
      if(!is.data.frame(esiAdducts)){
        stop('esiAdducts argument must be a character vector or data.frame.\n')
      }
      if(!all(c("name", "nmol", "Ch", "massDiff", "mode") %in% colnames(esiAdducts))){
        stop('The esiAdducts data.frame must contain the following column names:\n', paste0(1:5, c(". name", ". nmol", ". Ch", ". massDiff", ". mode"), collapse='\n'), '\n\nsee the adduct2mass function for further details.\n')
      }
    } else {
      if(Parameters(object)$mode == 'neg'){
        data(negESIAdducts)
        esiAdducts <- negESIAdducts
        message('default negative mode ', length(esiAdducts), ' adducts and in-source fragments.\n')
        flush.console()
      esiAdducts <- adduct2mass(esiAdducts)
      } else if(Parameters(object)$mode == 'pos'){
        data(posESIAdducts)
        esiAdducts <- posESIAdducts
        message('default positive mode ', length(esiAdducts), ' adducts and in-source fragments.\n')
        flush.console()
        esiAdducts <- adduct2mass(esiAdducts)
      }
    }
    indxTmp <- is.na(esiAdducts$massDiff)
    if(any(indxTmp)){
     warning('The adduct2mass function returned NAs for the following ESI adducts:\n',
             paste0(esiAdducts$name[indxTmp], collapse='\n'),
             '\nand will be removed.', immediate. = TRUE) 
      flush.console()
      esiAdducts <- esiAdducts[indxTmp == FALSE, , drop=F]
    }
    
    if(length(DBanno(object)) == 0){
     DBanno(object) <- vector('list', length(compSpectra(object)))
     names(DBanno(object)) <- names(compSpectra(object))
    }
    # feature subset
    if(is.null(featureSubSet)){
      featureSubSet <- names(compSpectra(object))
    }
    featureSubSet.indx <- featureSubSet %in% names(compSpectra(object)) 
    if(any(featureSubSet.indx == FALSE)){
      stop("The following composite spectra subset names do not match :", 
           paste0(sapply(featureSubSet[featureSubSet.indx == FALSE], message), "\n"))
    }
    # obtain MS1 mzs from compMS2 object
    unknowns <- sapply(metaData(object)[featureSubSet], function(x) unlist(x[grep("MS1_mz", names(x))])[1])
    if(!is.null(esiAdductsBySpec)){
    esiAdductsBySpec <- esiAdductsBySpec[featureSubSet]
    }
    # tmp substr mode indx
    if(!is.null(SubStrs)){
      if(SubStrs[1] == 'All'){
        SubStrs <- Substructure_masses
        mode.substr.indx <- SubStrs[, Parameters(object)$mode] == 1 & SubStrs$SubStructure == 1
        subStrMasses <- SubStrs$mass.shift[mode.substr.indx]
        subStrNames <- SubStrs$SubStructure_type[mode.substr.indx]
        # remove zero mass shifts
        mass.shift.zero <- subStrMasses != 0
        subStrMasses <- subStrMasses[mass.shift.zero]
        subStrNames <- subStrNames[mass.shift.zero]
        # duplicated substr types
        dupl.substr.type <- duplicated(subStrNames) == F
        subStrMasses <- subStrMasses[dupl.substr.type]
        subStrNames <- subStrNames[dupl.substr.type]
        names(subStrMasses) <- subStrNames  
      } else {
      
      stopifnot(is.numeric(SubStrs))
      if(is.null(names(SubStrs))){
        stop("Argument SubStrs is supplied but corresponding substructure names have not been set (i.e. SubStrs names attribute is null)")
      }
      subStrMasses <- SubStrs
     }
    } else {
    subStrMasses <- 0
    names(subStrMasses) <- ''
    }
     
    DBAnnoMatches <- monoMassMatch(unknowns = unknowns, metMasses.df =  metDB, 
                                   esiAdducts=esiAdducts,
                                   subStrMasses = subStrMasses, 
                                   mode = Parameters(object)$mode, 
                                   ppm = ppm,
                                   nCores = Parameters(object)$nCores)
    # subset by esi adduct type by spectrum
    if(!is.null(esiAdductsBySpec)){
      message('Filtering out possible false positives based on ESI adducts included in the MS1 feature table.\n')
      flush.console()
      befRem <- sum(sapply(DBAnnoMatches, nrow))
      for(j in 1:length(DBAnnoMatches)){
        DBannoTmp <- DBAnnoMatches[[j]]
          if(!is.null(DBannoTmp)){
            # if(esiAdductsBySpec[j] == ''){
            # esiDetTmp <- esiAdducts$name
            # } else {
            esiDetTmp <- strsplit(esiAdductsBySpec[j], ' ')[[1]]  
            # }
          esiDetTmp <- c(ifelse(Parameters(object)$mode == 'neg', '[M-H]-', '[M+H]+'), 
                         esiDetTmp)
          DBAnnoMatches[[j]] <- DBannoTmp[DBannoTmp$ESI_type %in% esiDetTmp, , drop=FALSE]
          }
      }
      aftRem <- sum(sapply(DBAnnoMatches, nrow))
      message(prettyNum(aftRem, big.mark = ','), ' annotations remaining following false positive removal.\n',
              round(100 - {{aftRem/befRem} * 100}, 2), '% of the annotations were removed.\n')
      flush.console()
    }
    # if neccessary and present include precalc MD and chemFP
    # if(includeMDChemFP == TRUE){
    # mdChemColIndx <- grepl('^MD_|chemFP', colnames(metDB))
    #   if(any(mdChemColIndx)){
    #   DBAnnoMatches <- lapply(DBAnnoMatches, function(x){
    #     indxTmp <- match(x$DBid, metDB$Unique_DB_ID) 
    #     return(cbind(x, metDB[indxTmp, mdChemColIndx, drop=FALSE]))
    #   }) 
    #   } 
    # }
    # any previous annotations
    #rbind previous DB search results
      indxTmp <- match(featureSubSet, names(DBanno(object)))
      DBAnnoMatches <- lapply(1:length(indxTmp), function(x){
        objIndx <- indxTmp[x]
        if(!is.null(DBanno(object)[[objIndx]])){
        # ensure all column names match else add empty NA 
        diffCols <- setdiff(colnames(DBanno(object)[[objIndx]]), colnames(DBAnnoMatches[[x]]))
        if(length(diffCols) > 0){
        emptyDf <- data.frame(matrix(NA, nrow=nrow(DBAnnoMatches[[x]]), 
                                     ncol=length(diffCols)), 
                              stringsAsFactors = FALSE)
        colnames(emptyDf) <- diffCols
        DBAnnoMatches[[x]] <- cbind(DBAnnoMatches[[x]], emptyDf)  
        }
        diffCols <- setdiff(colnames(DBAnnoMatches[[x]]), colnames(DBanno(object)[[objIndx]]))
        if(length(diffCols) > 0){
          emptyDf <- data.frame(matrix(NA, nrow=nrow(DBanno(object)[[objIndx]]), 
                                       ncol=length(diffCols)), 
                                stringsAsFactors = FALSE)
          colnames(emptyDf) <- diffCols
          DBanno(object)[[objIndx]] <- cbind(DBanno(object)[[objIndx]], emptyDf)  
        }
        tmp.df <- rbind(DBanno(object)[[objIndx]], DBAnnoMatches[[x]])
        # remove duplicates
        tmp.df <- tmp.df[duplicated(tmp.df$DBid) == FALSE, ]
        return(tmp.df)
      } else {
        return(DBAnnoMatches[[x]])
      }
      })
     
    names(DBAnnoMatches) <- featureSubSet
    DBanno(object)[featureSubSet] <- DBAnnoMatches
    return(object)
}) # end function
