#' Monoisotopic mass matching unknown to data base entry
#' @importFrom Rcpp evalCpp
#' @useDynLib compMS2Miner
monoMassMatch <- function(unknowns = NULL, metMasses.df = NULL, esiAdducts=NULL,
                          subStrMasses = NULL, mode="pos", ppm = 10, 
                          nCores=NULL){
  # error handling
  if(is.null(unknowns)){
    stop("unknowns argument is missing with no default")    
  } else if (is.null(esiAdducts)){
    stop("esiAdducts argument is missing with no default")    
  } else if (is.null(metMasses.df)){
    stop("metMasses.df argument is missing with no default")    
  } else {
    if(!all(c("Unique_DB_ID", "monoisotopic_weight", "name", "SMILES") %in% colnames(metMasses.df))){
      stop('metMasses.df colnames for unique database IDs, monoisotopic weights, names and canonical SMILES codes of data base
           entries must be named "Unique_DB_ID", "monoisotopic_weight", "name" and "SMILES" respectively' )
    }
    metMasses <- metMasses.df$monoisotopic_weight
    DB_unique_IDs <- metMasses.df$Unique_DB_ID
    DB_entry_names <- metMasses.df$name
    DB_SMILES <- metMasses.df$SMILES
    dbEntryWebAddress <- metMasses.df$WebAddress
    }  
  if(!grepl("pos|neg", mode)){
    stop('mode must be one of either "pos" or "neg"') 
  } else if (!is.null(subStrMasses)){
    stopifnot(is.numeric(subStrMasses))
    if(is.null(names(subStrMasses))){
      stop("Argument subStrMasses is supplied but corresponding substructure names have not been set (i.e. subStrMasses names attribute is null)")
     }
  }
  
  if(is.null(subStrMasses)){
    subStrMasses <- 0
    names(subStrMasses) <- ''
  } else {
    subStrMasses <- c(0, subStrMasses)
    names(subStrMasses)[1] <- ''
  }
  # remove duplicates
  subStrMasses <- subStrMasses[duplicated(names(subStrMasses)) == FALSE]
  esiAdducts$nmol <- as.numeric(esiAdducts$nmol)
  esiAdducts$Ch <- as.numeric(esiAdducts$Ch)
  esiAdducts$massDiff <- as.numeric(esiAdducts$massDiff)
    # foreach and dopar from foreach package
    Results <- vector('list', length(unknowns))
    pmt <- proc.time()  
    # single-threaded
    message("matching ", length(unknowns), " MS1 feature(s) to ", 
            length(subStrMasses), ' substructure type(s)...\n')
    flush.console()
  
    for(i in 1:length(subStrMasses)){
    subStrTypeTmp <- names(subStrMasses)[i]
    subStrMassTmp <- subStrMasses[i]
    ###calculate parent ions from metMasses supplied 
    # electrospray adducts taken from 
      metMassesTmp <- metMasses + subStrMassTmp
      metMassesAllEsi <- matrix(0, ncol=nrow(esiAdducts), nrow=length(metMassesTmp))
            
      for(m in 1:nrow(esiAdducts)){
      metMassesAllEsi[, m] <- {{esiAdducts$nmol[m] * metMassesTmp}/esiAdducts$Ch[m]}                               + esiAdducts$massDiff[m]
      }
      
      esiNamesTmp <- rep(esiAdducts$name, each=nrow(metMassesAllEsi))
      dbIdsTmp <- rep(1:nrow(metMasses.df), ncol(metMassesAllEsi))
      metMassesAllEsi <- as.vector(metMassesAllEsi)
     
        # max min range to limit search space
        uk.range <- range(unknowns)
        indxTmp <- metMassesAllEsi >= (uk.range[1]-1) & metMassesAllEsi <= (uk.range[2]+1)
          metMassesAllEsiSub <- metMassesAllEsi[indxTmp]
          esiNamesTmpSub <- esiNamesTmp[indxTmp]
          dbIdsTmpSub <- dbIdsTmp[indxTmp] 
          message("matching ", length(unknowns), " unknown(s) to ", 
                  format(length(metMassesAllEsiSub), big.mark=",", 
                         scientific=FALSE),
                  " possible ESI artefact/ in-source fragment mass shifts ", ifelse(subStrTypeTmp == '', '(no substructure)', paste0('and the following substructure mass shift: ', subStrTypeTmp)), "...please wait\n")
          flush.console()
          
          matches.tmp <- vector("list", length(unknowns))
          matches.tmp <- compMS2Miner:::exactMassMatchCpp(matches.tmp, unknowns, metMassesAllEsiSub, 1:length(metMassesAllEsiSub), ppm, 1E06)
          
          matches.tmp <- lapply(1:length(matches.tmp), function(x){
            mTmp <- as.numeric(matches.tmp[[x]])
            monoMass <- metMassesAllEsiSub[mTmp]
            ESI_type <- esiNamesTmpSub[mTmp]
            SubStr_type <- rep(subStrTypeTmp, length(mTmp))
            indx.tmp <- dbIdsTmpSub[mTmp]
            tmp.df <- data.frame(WebAddress = dbEntryWebAddress[indx.tmp],
                                 DBid = DB_unique_IDs[indx.tmp], 
                                 DBname = DB_entry_names[indx.tmp],
                                 SMILES = DB_SMILES[indx.tmp],
                                 ESI_type = gsub(';', '', ESI_type), 
                                 SubStr_type = SubStr_type,
                                 observedMass =  rep(unknowns[x], length(monoMass)),
                                 candidateMass = as.numeric(monoMass), stringsAsFactors = FALSE)
            tmp.df$ppmMatch <- ((tmp.df$observedMass - tmp.df$candidateMass) / tmp.df$observedMass) * 1E06
            tmp.df$BestAnno <- rep(FALSE, nrow(tmp.df))
            return(tmp.df)})        
          
        Results <- lapply(1:length(Results), function(x){
                          rbind(Results[[x]], matches.tmp[[x]])})
        }
        proc.time() - pmt
        # 320 secs 198 unknowns 18 substructure mass shifts.
        # nMatches 
        nMatches <- sapply(Results, nrow)
        message(format(sum(nMatches), big.mark=",", scientific=FALSE), ' database matches were made to ', length(Results), ' unknown MS1 precursor(s) (mean=', format(round(sum(nMatches)/length(Results), 0), big.mark=",", scientific=FALSE), ', min=', format(min(nMatches), big.mark=",", scientific=FALSE), ', max=', format(max(nMatches), big.mark=",", scientific=FALSE), ').\n', nrow(esiAdducts), ' different ESI adducts/ in-source fragments', ifelse(all(names(subStrMasses) %in% ''), ' were considered...\n', paste0(' + ', length(subStrMasses), ' substructure mass shifts were considered...\n')))
        flush.console()
        
        return(Results)
        # TO DO: inverse mass shift calc.
} # end function
