load(file='compMS2object.RData')
# if necessary read in r script and sessioninfo text files
# suppressWarnings(suppressMessages(library(compMS2Miner)))
suppressWarnings(suppressMessages(library(igraph)))
suppressWarnings(suppressMessages(library(rhandsontable)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(shiny)))
# read session info if present
seshInfo <- list.files(pattern='sessionInfo_')
if(length(seshInfo) > 0){
  if(length(seshInfo) > 1){
    datesSum <- gsub('sessionInfo_|\\.txt$', '', basename(seshInfo))
    datesSum <- sapply(strsplit(datesSum, ''), function(x) sum(as.numeric(x)))
    seshInfo <- seshInfo[which.max(datesSum)]
  }
seshTxtTmp <- c(paste0('<h5>', seshInfo, '</h5>', '</br>'), 
                            paste(readLines(seshInfo), '</br>'))  
} else {
  # generate sessionInfo text file
  seshInfoTxt <- paste0('sessionInfo_', gsub('-', '', Sys.Date()),
                        '.txt')
  writeLines(capture.output(sessionInfo()), seshInfoTxt)
  seshInfo <- list.files(pattern='sessionInfo_')
  seshTxtTmp <- c(paste0('<h5>', seshInfo, '</h5>', '</br>'), 
                  paste(readLines(seshInfo), '</br>'))  
} 
# see if there is a pdf to display
pdfFile <- list.files(path='www', pattern='\\.pdf$')[1]
if(length(pdfFile) == 0){
  # pdfFile <- 'https://github.com/WMBEdmands/compMS2Miner/blob/master/vignettes/CompMS2miner_Workflow.pdf' 
  pdfFile <- 'https://raw.githubusercontent.com/WMBEdmands/compMS2Miner/master/vignettes/CompMS2miner_Workflow.pdf'
  # not found works just as well
}

composite_spectra <- compSpectra(object)
Features.v <- names(composite_spectra)
###DB search names
tmp.DBanno.res <- DBanno(object)

### best anno
tmp.BestAnno <- BestAnno(object)
metaData.tmp <- metaData(object)
UserComments.v <- vector("list", length(composite_spectra))
# best substructure anno
subStrAnno.df <- subStrAnno(object)

# metFrag 
tmp.metFrag <- MetFrag(object)
# CFM 
tmp.CFM <- CFM(object)
###order by EIC number
EICorderIndx <- order(as.numeric(gsub(".+_","",Features.v)))
Features.v <- Features.v[EICorderIndx]
composite_spectra <- composite_spectra[EICorderIndx]
metaData.tmp <- metaData.tmp[EICorderIndx]

names(UserComments.v) <- names(composite_spectra)
if(length(tmp.DBanno.res) > 0){
  tmp.DBanno.res <- tmp.DBanno.res[EICorderIndx]
}
if(length(tmp.BestAnno) > 0){
  tmp.BestAnno <- tmp.BestAnno[EICorderIndx]
}
if(length(tmp.metFrag) > 0){
  tmp.metFrag <- tmp.metFrag[EICorderIndx]
}
if(length(tmp.CFM) > 0){
  tmp.CFM <- tmp.CFM[EICorderIndx]
  cfmSelectTable <- sapply(tmp.CFM, function(x){
     if(!is.null(x)){
     indxDBid <- duplicated(x$DBid) == F
     nPeaksEx <- table(x$DBid)
     x <- x[indxDBid, c('DBid', 'WebAddress', 'DBname', 'CFM_totPropEx'), drop=F]
     x$nPeaksEx[match(x$DBid, names(nPeaksEx))] <- nPeaksEx
     return(x)
     } else {
       return(NULL)
     }})
} else {
  cfmSelectTable <- NULL
}
# DB match df to string match
if(length(DBanno(object)) > 0){
  DBmatches <-  t(sapply(tmp.DBanno.res, function(x){
    if(!is.null(x)){
    names.tmp <- x$DBname
    ESI_type.tmp <- x$ESI_type
    SubStr_type.tmp <- x$SubStr_type
    return(data.frame(names.tmp, ESI_type.tmp, SubStr_type.tmp, 
                            stringsAsFactors = F))
    } else {
      return(data.frame(names.tmp='', ESI_type.tmp='', SubStr_type.tmp='', 
                 stringsAsFactors = F))
    }
  }))
} else {
  DBmatches <- matrix("The metID.dbAnnotate function has not yet been run", 
                      ncol = 2, nrow = length(composite_spectra))
}

# DB match df to string match
if(length(BestAnno(object)) > 0){
  DBBestMatches <-  t(sapply(tmp.BestAnno, function(x){
    if(!is.null(x)){
      names.tmp <- x$DBname
      ESI_type.tmp <- x$ESI_type
      SubStr_type.tmp <- x$SubStr_type
      return(data.frame(names.tmp, ESI_type.tmp, SubStr_type.tmp, 
                        stringsAsFactors = F))
    } else {
      return(data.frame(names.tmp='', ESI_type.tmp='', SubStr_type.tmp='', 
                        stringsAsFactors = F))
    }
  }))
} else {
  DBBestMatches <- matrix("The metID.dbProb function has not yet been run", 
                          ncol = 2, nrow = length(composite_spectra))
}
# substructures id'd
substrIded <- sapply(composite_spectra, function(x){
  if(is.data.frame(x)){
    'Frag.ID' %in% colnames(x)
  } else {
  FALSE  
  }})
if(any(substrIded)){
  
  SubStr_types <- t(sapply(composite_spectra, function(x){
    if(!is.null(x)){
      Frag.ID.tmp <- paste0("Frag_", as.character(x$Frag.ID))
      Neutral.loss.tmp <- paste0("NeutLoss_", as.character(x$Neutral.loss))
      SubStrType.tmp <- data.frame(Frag.ID.tmp, Neutral.loss.tmp, 
                                   stringsAsFactors = F)
    } else {
      return(data.frame(Frag.ID.tmp='', Neutral.loss.tmp='', 
                        stringsAsFactors = F))
    }
  }))
} else {
  SubStr_types <- matrix("The subStructure.Annotate function has not yet been run", 
                         ncol = 2, nrow = length(composite_spectra))
}
###obtain SubStr types for filtration
SubStrType.inputs <- c(unique(unlist(SubStr_types[, 1])), 
                       unique(unlist(SubStr_types[, 2])))

SubStrType.inputs <- SubStrType.inputs[-grep("Frag_$|NeutLoss_$", SubStrType.inputs)]

if(all(SubStrType.inputs == "The subStructure.Annotate function has not yet been run")){
  SubStrType.inputs <- "The subStructure.Annotate function has not yet been run"
}
###extract mass and RT values
mass.v <- sapply(metaData.tmp,function(x) as.numeric(x[grep("MS1_mz", names(x))][[1]][1]))
RT.v <- sapply(metaData.tmp, function(x) as.numeric(x[grep("MS1_RT", names(x))][[1]][1]))
meanPrecursorInt <- sapply(metaData.tmp, function(x) mean(as.numeric(unlist(x[grep("precursorIntensity", names(x))]))))
allFeatTable <- data.frame(specNames=names(mass.v), mass=mass.v, rt=RT.v, meanPrecursorInt=meanPrecursorInt,  stringsAsFactors = F)
allFeatTable$meanPrecursorInt[is.na(allFeatTable$meanPrecursorInt)] <- 1
if(all(meanPrecursorInt == 0)){
allFeatTable$precursorInt_group <- 1  
} else {
allFeatTable$precursorInt_group <- as.numeric(cut(meanPrecursorInt, 10))
}
#  

# substr annotations if necessary
if(nrow(subStrAnno.df) > 0){
  # sort by EICno
  subStrAnno.df <- subStrAnno.df[order(as.numeric(gsub(".+_","", subStrAnno.df$compSpecName))), , drop=F]
  duplEntTmp <- duplicated(subStrAnno.df$compSpecName) == F
  bestAnnoSubStr <- subStrAnno.df[duplEntTmp, , drop=F]
  matchIndxTmp <- match(allFeatTable$specNames, bestAnnoSubStr$compSpecName)
  bestAnnoSubStr <- bestAnnoSubStr[matchIndxTmp, c(1, 3:5)]
  colnames(bestAnnoSubStr)[1] <- 'possible substructure'
  bestAnnoSubStr$SumRelInt <- suppressWarnings(ifelse(is.na(bestAnnoSubStr$SumRelInt), NA, round(as.numeric(bestAnnoSubStr$SumRelInt), 2))) 
  allFeatTable <- cbind(allFeatTable, bestAnnoSubStr)
  subStrAnno.df$SumRelInt <- ifelse(subStrAnno.df$SumRelInt == "no substructure detected", 
                                    0, subStrAnno.df$SumRelInt)
  subStrAnno.df$SumRelInt <- round(as.numeric(subStrAnno.df$SumRelInt), digits=2)  
}
subStrAnno.list <- split(subStrAnno.df, subStrAnno.df$compSpecName)
# sort 
subStrAnno.list <- subStrAnno.list[order(as.numeric(gsub(".+_","", names(subStrAnno.list))))]
subStrAnno.inputs <- unique(subStrAnno.df$SubStrType)

TotalFeatures <- length(unique(gsub(".+_", "", Features.v)))
TotalCompSpectra <- length(Features.v)
# network graph
if(!is.null(network(object)$corrNetworkGraph)){
corrNetTmp <- network(object)$corrNetworkGraph
corrLayoutTmp <- network(object)$corrLayout
colnames(corrLayoutTmp)[1:2] <- c('xvar', 'yvar')
# rescale layout
corrScaledLayout <- data.frame(corrLayoutTmp, stringsAsFactors = F)
corrScaledLayout$xvar <- scales::rescale(corrScaledLayout$xvar, c(-1, 1))
corrScaledLayout$yvar <- scales::rescale(corrScaledLayout$yvar, c(-1, 1))
corrScaledLayout[, 4] <- round(corrScaledLayout[, 4], 4)
corrScaledLayout[, 5] <- round(corrScaledLayout[, 5], 4)
corrNetMatchIndx <- match(igraph::V(corrNetTmp)$name, Features.v)
# add possible substructure column
possible_substructures <- ifelse(grepl('^CC_', Features.v[corrNetMatchIndx]), allFeatTable[corrNetMatchIndx[!is.na(corrNetMatchIndx)], 'possible substructure'], 'no composite spectrum')
corrScaledLayout <- cbind(corrScaledLayout, possible_substructures)

igraph::V(corrNetTmp)$MS2netColours <- ifelse(grepl('^CC_', Features.v[corrNetMatchIndx]), "#D55E00", "#0072B2")
igraph::V(corrNetTmp)$vertexShapes <- rep('circle', length(corrNetMatchIndx))
igraph::V(corrNetTmp)$vertexSize <- rep(4, length(corrNetMatchIndx))

}
# spectral similarity
if(!is.null(network(object)$specSimGraph)){
  specSimNetTmp <- network(object)$specSimGraph
  specSimLayoutTmp <- network(object)$specSimLayout
  colnames(specSimLayoutTmp)[1:2] <- c('xvar', 'yvar')
  # rescale layout
 specSimScaledLayout <- data.frame(specSimLayoutTmp, stringsAsFactors = F)
 specSimScaledLayout$xvar <- scales::rescale(specSimScaledLayout$xvar, c(-1, 1))
 specSimScaledLayout$yvar <- scales::rescale(specSimScaledLayout$yvar, c(-1, 1))
 specSimScaledLayout$name <- igraph::V(specSimNetTmp)$name
 scIdx <- match(igraph::V(specSimNetTmp)$name, Features.v)
 
 specSimScaledLayout$mass <- 0
 specSimScaledLayout$rt <- 0
 specSimScaledLayout$mass[!is.na(scIdx)] <- round(allFeatTable$mass[scIdx[!is.na(scIdx)]], 4)
 specSimScaledLayout$rt[!is.na(scIdx)] <- round(allFeatTable$rt[scIdx[!is.na(scIdx)]], 2)
 igraph::V(specSimNetTmp)$MS2netColours <- ifelse(!is.na(scIdx), "#D55E00", "#0072B2")
 igraph::V(specSimNetTmp)$vertexShapes <- rep('circle', length(scIdx))
 igraph::V(specSimNetTmp)$vertexSize <- rep(4, length(scIdx))
}

# reconsubstructure similarity
if(!is.null(network(object)$reconSubGraph)){
  reconSubNetTmp <- network(object)$reconSubGraph
  reconSubLayoutTmp <- network(object)$reconSubLayout
  colnames(reconSubLayoutTmp)[1:2] <- c('xvar', 'yvar')
  # rescale layout
  reconSubScaledLayout <- data.frame(reconSubLayoutTmp, stringsAsFactors = F)
  reconSubScaledLayout$xvar <- scales::rescale(reconSubScaledLayout$xvar, c(-1, 1))
  reconSubScaledLayout$yvar <- scales::rescale(reconSubScaledLayout$yvar, c(-1, 1))
  reconSubScaledLayout$name <- igraph::V(reconSubNetTmp)$name
  rcIdx <- match(igraph::V(reconSubNetTmp)$name, Features.v)
  reconSubScaledLayout$mass <- 0
  reconSubScaledLayout$rt <- 0
  reconSubScaledLayout$mass[!is.na(rcIdx)] <- round(allFeatTable$mass[rcIdx[!is.na(rcIdx)]], 4)
  reconSubScaledLayout$rt[!is.na(rcIdx)] <- round(allFeatTable$rt[rcIdx[!is.na(rcIdx)]], 2)
  igraph::V(reconSubNetTmp)$vertexShapes <- rep('circle', length(rcIdx))
}
# chemical similarity
# if(!is.null(network(object)$chemSimGraph)){
#   specSimNetTmp <- network(object)$specSimGraph
#   specSimLayoutTmp <- network(object)$specSimLayout
#   colnames(specSimLayoutTmp)[1:2] <- c('xvar', 'yvar')
#   # rescale layout
#   specSimScaledLayout <- data.frame(specSimLayoutTmp, stringsAsFactors = F)
#   specSimScaledLayout$xvar <- scales::rescale(specSimScaledLayout$xvar, c(-1, 1))
#   specSimScaledLayout$yvar <- scales::rescale(specSimScaledLayout$yvar, c(-1, 1))
#   specSimScaledLayout$name <- igraph::V(specSimNetTmp)$name
#   # specSimScaledLayout$mzmed <- round(specSimScaledLayout$mzmed, 4)
#   # specSimScaledLayout$rtmed <- round(specSimScaledLayout$rtmed, 4)
#   specSimMatchIndx <- match(igraph::V(specSimNetTmp)$name, Features.v)
#   specSimScaledLayout <- cbind(specSimScaledLayout, allFeatTable[specSimMatchIndx, c('mass', 'rt', 'possible substructure')]) 
#   igraph::V(specSimNetTmp)$MS2netColours <- ifelse(!is.na(specSimMatchIndx), "#D55E00", "#0072B2")
#   igraph::V(specSimNetTmp)$vertexShapes <- rep('circle', length(specSimMatchIndx))
#   igraph::V(specSimNetTmp)$vertexSize <- rep(4, length(specSimMatchIndx))
# }
# 
# # if all 3 graphs then create combined
# allGraphsIndx <- c('corrNetworkGraph', 'specSimGraph', 'chemSimGraph') %in% names(network(object))
# if(all()){
#   
# }

# met Id comments table
metIDcomments <- Comments(object)


# spectral database indx
if(length(object@spectralDB) > 0){
indxSpectralDb <- sapply(object@spectralDB, length) > 0
indxSpectralDb <- which(indxSpectralDb[EICorderIndx])
specDBmatches <- object@spectralDB[EICorderIndx]
} else {
indxSpectralDb <- as.numeric() 
}
# in silico fragment indx
if(length(object@inSilico) > 0){
  inSilicoIndx <- rep(FALSE, length(Features.v))
  if(!is.null(MetFrag(object))){
  inSilicoIndx[sapply(MetFrag(object), function(x){
     if(!is.null(x)){ifelse(ncol(x) > 1, TRUE, FALSE)} else {F}})] <- T
  }
  if(!is.null(CFM(object))){
    inSilicoIndx[sapply(CFM(object), function(x){
      if(!is.null(x)){ifelse(ncol(x) > 1, TRUE, FALSE)} else {F}})] <- T
  }
  inSilicoIndx <- which(inSilicoIndx[EICorderIndx])
} else {
  inSilicoIndx <- as.numeric() 
}
# retention time prediction model
rtPredIndx <- sapply(tmp.BestAnno, function(x) any(grepl('predRtDev', colnames(x))))

if(any(rtPredIndx)){
  specNamesTmp <- unlist(mapply(rep, names(tmp.BestAnno), each=sapply(tmp.BestAnno, nrow))) 
  bestAnnoDf <- do.call(rbind, tmp.BestAnno)
  bestAnnoDf <- bestAnnoDf[, c("WebAddress", "DBid", "DBname", "SMILES",
                               "SubStr_type", "observedMass", "trainingSet", "predRts", "predRtDev", 'ms1Rt'), drop=F]
  bestAnnoDf <- cbind(specNames=specNamesTmp, bestAnnoDf)
  bestAnnoDf <- bestAnnoDf[bestAnnoDf$predRtDev != '', ]
  
  if(!is.null(object@rtPred$standardsTable)){
  stdTableTmp <- object@rtPred$standardsTable
  colnames(stdTableTmp) <- c('DBname', 'SMILES', 'ms1Rt', 'predRts')
  stdTableTmp$predRtDev <- stdTableTmp$predRts - stdTableTmp$ms1Rt
  stdTableTmp$trainingSet <- TRUE
  stdTableTmp$WebAddress <- ''
  stdTableTmp$DBid <- ''
  stdTableTmp$SubStr_type <- 'standard'
  stdTableTmp$observedMass <- ''
  stdTableTmp$specNames <- 'standard'
  stdTableTmp <- stdTableTmp[, match(colnames(bestAnnoDf), colnames(stdTableTmp))]
  bestAnnoDf <- rbind(bestAnnoDf, stdTableTmp)
  }
  indxTmp <- bestAnnoDf$WebAddress != ''
  bestAnnoDf$DBid[indxTmp]  <- paste0("<a href='http://", bestAnnoDf$WebAddress[indxTmp],  
                            bestAnnoDf$DBid[indxTmp], "' target='_blank'>",  
                            bestAnnoDf$DBid[indxTmp],  "</a>")
  bestAnnoDf$SMILES <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", 
                             gsub("/|\\\\", "",  bestAnnoDf$SMILES), "/PNG", 
                             "' target='_blank'>", paste0(substring(bestAnnoDf$SMILES,  1,  6),  "..."),
                             "</a>")
 bestAnnoDf$WebAddress <- NULL
 bestAnnoDf$SubStr_type <- NULL
 bestAnnoDf$predRtDev <- round(as.numeric(bestAnnoDf$predRtDev), 3)
 bestAnnoDf$predRts <- round(as.numeric(bestAnnoDf$predRts), 3)
 bestAnnoDf$ms1Rt <- as.numeric(bestAnnoDf$ms1Rt)
 lmPredRt <- lm(bestAnnoDf$predRts[as.logical(bestAnnoDf$trainingSet)] ~ bestAnnoDf$ms1Rt[as.logical(bestAnnoDf$trainingSet)])
 }
