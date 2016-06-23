load(file='compMS2object.RData')

suppressWarnings(suppressMessages(library(CompMS2miner)))
suppressWarnings(suppressMessages(library(igraph)))
suppressWarnings(suppressMessages(library(rhandsontable)))
suppressWarnings(suppressMessages(library(scales)))

composite_spectra <- object@compSpectra
Features.v <- names(composite_spectra)
###DB search names
tmp.DBanno.res <- object@DBanno

### best anno
tmp.BestAnno <- object@BestAnno
metaData.tmp <- object@metaData
UserComments.v <- vector("list", length(composite_spectra))
# best substructure anno
subStrAnno.df <- object@subStrAnno

# metFrag 
tmp.metFrag <- object@MetFrag
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
  tmp.metFrag<- tmp.metFrag[EICorderIndx]
}
# DB match df to string match
if(length(object@DBanno) > 0){
  DBmatches <-  t(sapply(tmp.DBanno.res, function(x){
    names.tmp <- x$DBname
    ESI_type.tmp <- x$ESI_type
    SubStr_type.tmp <- x$SubStr_type
    names.tmp <- data.frame(names.tmp, ESI_type.tmp, SubStr_type.tmp, 
                            stringsAsFactors = F)
  }))
} else {
  DBmatches <- matrix("The metID.dbAnnotate function has not yet been run", 
                      ncol = 2, nrow = length(composite_spectra))
}

# DB match df to string match
if(length(object@BestAnno) > 0){
  DBBestMatches <-  t(sapply(tmp.BestAnno, function(x){
    names.tmp <- x$DBname
    ESI_type.tmp <- x$ESI_type
    SubStr_type.tmp <- x$SubStr_type
    names.tmp <- data.frame(names.tmp, ESI_type.tmp, SubStr_type.tmp, 
                            stringsAsFactors = F)
  }))
} else {
  DBBestMatches <- matrix("The metID.dbProb function has not yet been run", 
                          ncol = 2, nrow = length(composite_spectra))
}
# substructures id'd
if(!is.null(composite_spectra[[1]]$Frag.ID)){
  
  SubStr_types <- t(sapply(composite_spectra, function(x){
    Frag.ID.tmp <- paste0("Frag_", as.character(x$Frag.ID))
    Neutral.loss.tmp <- paste0("NeutLoss_", as.character(x$Neutral.loss))
    SubStrType.tmp <- data.frame(Frag.ID.tmp, Neutral.loss.tmp, 
                                 stringsAsFactors = F)
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
allFeatTable <- data.frame(specNames=names(mass.v), mass=mass.v, rt=RT.v, meanPrecursorInt=meanPrecursorInt, precursorInt_group=as.numeric(cut(meanPrecursorInt, 10)),  stringsAsFactors = F) #  

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
if(!is.null(object@network$corrNetworkGraph)){
corrNetTmp <- object@network$corrNetworkGraph
corrLayoutTmp <- object@network$corrLayout
colnames(corrLayoutTmp)[1:2] <- c('xvar', 'yvar')
# rescale layout
corrScaledLayout <- data.frame(corrLayoutTmp, stringsAsFactors = F)
corrScaledLayout$xvar <- scales::rescale(corrScaledLayout$xvar, c(-1, 1))
corrScaledLayout$yvar <- scales::rescale(corrScaledLayout$yvar, c(-1, 1))
corrScaledLayout$mzmed <- round(corrScaledLayout$mzmed, 4)
corrScaledLayout$rtmed <- round(corrScaledLayout$rtmed, 4)
corrNetMatchIndx <- match(as.numeric(igraph::V(corrNetTmp)$name), as.numeric(gsub(".+_", "", Features.v)))
# add possible substructure column
possible_substructures <- ifelse(!is.na(corrNetMatchIndx), allFeatTable[corrNetMatchIndx[!is.na(corrNetMatchIndx)], 'possible substructure'], 'no composite spectrum')
corrScaledLayout <- cbind(corrScaledLayout, possible_substructures)

igraph::V(corrNetTmp)$MS2netColours <- ifelse(!is.na(corrNetMatchIndx), "#D55E00", "#0072B2")
igraph::V(corrNetTmp)$vertexShapes <- rep('circle', length(corrNetMatchIndx))
igraph::V(corrNetTmp)$vertexSize <- rep(4, length(corrNetMatchIndx))

}
# spectral similarity
if(!is.null(object@network$specSimGraph)){
  specSimNetTmp <- object@network$specSimGraph
  specSimLayoutTmp <- object@network$specSimLayout
  colnames(specSimLayoutTmp)[1:2] <- c('xvar', 'yvar')
  # rescale layout
 specSimScaledLayout <- data.frame(specSimLayoutTmp, stringsAsFactors = F)
 specSimScaledLayout$xvar <- scales::rescale(specSimScaledLayout$xvar, c(-1, 1))
 specSimScaledLayout$yvar <- scales::rescale(specSimScaledLayout$yvar, c(-1, 1))
 specSimScaledLayout$name <- igraph::V(specSimNetTmp)$name
 # specSimScaledLayout$mzmed <- round(specSimScaledLayout$mzmed, 4)
 # specSimScaledLayout$rtmed <- round(specSimScaledLayout$rtmed, 4)
  specSimMatchIndx <- match(igraph::V(specSimNetTmp)$name, Features.v)
  specSimScaledLayout <- cbind(specSimScaledLayout, allFeatTable[specSimMatchIndx, c('mass', 'rt', 'possible substructure')]) 
  igraph::V(specSimNetTmp)$MS2netColours <- ifelse(!is.na(specSimMatchIndx), "#D55E00", "#0072B2")
  igraph::V(specSimNetTmp)$vertexShapes <- rep('circle', length(specSimMatchIndx))
  igraph::V(specSimNetTmp)$vertexSize <- rep(4, length(specSimMatchIndx))
}

# met Id comments table
if(nrow(object@Comments) > 0){
metIDcomments <- object@Comments
} else {
metIDcomments <- data.frame(compSpectrum=Features.v, possible_identity=rep('', length(Features.v)), compound_class=rep('', length(Features.v)), user_comments=rep('', length(Features.v)), stringsAsFactors = F)
}

# spectral database indx
if(length(object@spectralDB) > 0){
indxSpectralDb <- sapply(object@spectralDB, length) > 0
indxSpectralDb <- which(indxSpectralDb[EICorderIndx])
specDBmatches <- object@spectralDB[EICorderIndx]
} else {
indxSpectralDb <- as.numeric() 
}