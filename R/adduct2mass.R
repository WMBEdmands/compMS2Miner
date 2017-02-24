#' converts ESI adduct names into table of monoisotopic masses
#' 
#' @details this function can be used to generate a table of ESI adducts for the
#' \code{\link{metID.dbAnnotate}} function.
#' @param adductNames character vector of ESI adduct names in a specific form for example see default names c('[M-H]-', '[2M+2CH3OH]2-', '[M-H+C2H4O2+Na]-'). The 
#' function was developed and tested against the 128 and 133 different ESI 
#' adducts and in-source fragment names for accuracy.
#' @return data.frame of n rows corresponding to each adduct name and 5 columns:
#' \enumerate{
#' \item "name" adduct name    
#' \item "nmol" number of molecules (e.g. 2M-H = 2, M-H = 1)   
#' \item "Ch" charge state      
#' \item "massDiff" summed monoisotopic mass difference.
#' \item "mode" ion polarity based on adduct name (e.g. "]-" = 'neg", and "]+" = 'pos")
#' }
#' @source element monoisotopic masses and natural abundances taken from \url{http://www.sisweb.com/referenc/source/exactmas.htm}, see ?exactMassEle.
#' @references 
#' \enumerate{
#' \item Stanstrup, J., Gerlich, M., Dragsted, L.O. et al. 
#' Anal Bioanal Chem (2013) 405: 5037. doi:10.1007/s00216-013-6954-6
#' }
#' @export
adduct2mass <- function(adductNames=c('[M-H]-', '[2M+2CH3OH]2-', '[M-H+C2H4O2+Na]-')){
  # create empty formula from exactMassEle
  # emptyForm <- as.list(rep(0, nrow(exactMassEle)))
  # names(emptyForm) <- exactMassEle$eleSymbol
  adductNames <- ifelse(adductNames == '[M+Hac-H]-', "[M-H+CH3COOH]-", adductNames)
  adductNames <- ifelse(adductNames == '[M+FA-H]-', "'[M-H+HCOOH]-'", adductNames)
  
  adductNamesTmp <- adductNames
  massElectron <- 0.00054857990924
  # remove comment strs e.g.  (McLafferty)
  adductNames <- gsub('[[:blank:]]\\(.+\\)', '', adductNames)
  # if charged and positive mode
  chIndx <- grepl('[0-9]\\+\\]', adductNames)
  if(any(chIndx)){
    nElectrons <- gsub('\\+\\].+', '', adductNames[chIndx])
    nElectrons <- as.numeric(gsub('.+[^0-9]', '', nElectrons))
    adductNames[chIndx] <- gsub('[0-9]+\\+\\]', ']', adductNames[chIndx])
    chIndx[grepl('\\-$', adductNames)] <- FALSE
  }
  
  resTmp <- lapply(strsplit(adductNames, '\\+|\\-|\\[|\\]'), function(x){
        nmol <- ifelse(x[2] == "M", 1, as.numeric(gsub('M', '', x[2])))
        ch <- ifelse(x[length(x)] == '', 1, as.numeric(x[length(x)]))
        x <- x[c(-1, -2, -length(x))]
        # if in brackets
        brIndx <- grepl('\\(', x)
        if(any(brIndx)){
        nCtmp <- regexpr('\\)', x[brIndx]) 
        nCtmp <- as.numeric(substr(x[brIndx], nCtmp + 1, nCtmp + 1))
        tmpRepBr <- paste0(rep(gsub('\\(|\\)[0-9].+', '', x[brIndx]), nCtmp), collapse = '') 
        x[brIndx] <- paste0(tmpRepBr, gsub('\\(.+\\)[0-9]', '', x[brIndx]), collapse = '') 
        }
        multAdd <- substring(x, 1, 1)
        multAddIndx <- grepl('[0-9]', multAdd)
        x <- ifelse(multAddIndx, substring(x, 2, nchar(x)), x)
        x <- x[x != '']
        adductMonoMassesTmp <- sapply(strsplit(x, ''), function(y){
          lowCaseTmp <- grep('[[:lower:]]', y)
          if(length(lowCaseTmp) > 0){
          y[lowCaseTmp - 1] <- paste0(y[lowCaseTmp - 1], y[lowCaseTmp]) 
          y <- y[-lowCaseTmp]
          }
          # n ele
          numTmp <- grep('[0-9]', y)
          if(length(numTmp) > 0){
          diffNum <- which(diff(numTmp) == 1) 
          if(length(diffNum) > 0){
          y[numTmp][diffNum] <- paste0(y[numTmp][diffNum], y[numTmp][diffNum + 1], collapse = '') 
          y <- y[-numTmp[diffNum + 1]]
          numTmp <- numTmp[-(diffNum + 1)]
          }
          y <- c(y, rep(y[numTmp - 1], as.numeric(y[numTmp]) - 1))
          y <- y[-numTmp]
          }
          tabY <- table(y)
          indxTmp <- match(exactMassEle$eleSymbol, names(tabY))
          # calculate summed monoiso mass based on natural abundance
          massAddTmp <- sum(sapply(which(!is.na(indxTmp)), function(z){
            monoMassesTmp <- as.numeric(strsplit(exactMassEle$monoMass[z], ' ')[[1]])
            natAbundTmp <- as.numeric(strsplit(exactMassEle$natAbund[z], ' ')[[1]])
            monoMassesTmp <- monoMassesTmp[which.max(natAbundTmp)]
            monoMassesTmp <- monoMassesTmp * tabY[indxTmp[z]]
            # natAbundTmp <- as.numeric(strsplit(exactMassEle$natAbund[z], ' ')[[1]])
            # if(length(monoMassesTmp) > 1){
            #    return(sum(sample(monoMassesTmp, size=tabY[indxTmp[z]], 
            #                      replace=TRUE, prob=natAbundTmp)))
            #   } else {
              return(monoMassesTmp)       
            # }
            }))
       }) # end mass calc
        if(any(multAddIndx)){
        adductMonoMassesTmp[multAddIndx] <- adductMonoMassesTmp[multAddIndx] * as.numeric(multAdd[multAddIndx])  
        }
        return(c(nmol, ch, paste0(adductMonoMassesTmp, collapse=' ')))
      })
   
  resTmpDf <- data.frame(do.call(rbind, resTmp), stringsAsFactors=FALSE)
  colnames(resTmpDf) <- c('nmol', 'Ch', 'massDiff')
  plusMinTmp <- strsplit(adductNames, '\\[|\\]|[A-Z]|[a-z]|[0-9]|\\)|\\(')
  for(j in 1:nrow(resTmpDf)){
  nElec <-  massElectron * as.numeric(resTmpDf$Ch[j])
  tmpPM <- plusMinTmp[[j]]
  tmpPM <- tmpPM[tmpPM != '']
  resTmpDf$mode[j] <- ifelse(tmpPM[length(tmpPM)] == "+", 'pos', 'neg')
  tmpPM <- tmpPM[-length(tmpPM)]
  if(length(tmpPM) > 1){
  # massesV <- c(0, as.numeric(strsplit(resTmpDf$massDiff[j], ' ')[[1]]))
  massesV <- as.numeric(strsplit(resTmpDf$massDiff[j], ' ')[[1]])
  subtrMasses <- sum(massesV[tmpPM == '-'], ifelse(resTmpDf$mode[j] == 'neg', -nElec, 0))
  addiMasses <- sum(massesV[tmpPM == '+'], ifelse(resTmpDf$mode[j] == 'pos', -nElec, 0), ifelse(chIndx[j], nElectrons[j] * massElectron, 0))
 formTmp <- paste0(0, ifelse(length(subtrMasses) > 0, paste0("-", subtrMasses), ''),
                   ifelse(length(addiMasses) > 0, paste0("+", addiMasses), ''),
                   collapse='')
  # add the subtraction of the correct number of electrons
  massDiffTmp <- eval(parse(text=formTmp)) 
  } else {
  formTmp <- paste0(0, tmpPM[1], sum(as.numeric(resTmpDf$massDiff[j]), 
                    ifelse(chIndx[j], nElectrons[j] * massElectron, 0)))
  massDiffTmp <- eval(parse(text=formTmp))
  massDiffTmp <- ifelse(massDiffTmp < 0, massDiffTmp + nElec, massDiffTmp - nElec)
  }
  resTmpDf$massDiff[j] <- massDiffTmp
  }
  resTmpDf <- cbind(name=adductNames, resTmpDf)
  resTmpDf$name <- adductNamesTmp
  return(resTmpDf)
} # end function
