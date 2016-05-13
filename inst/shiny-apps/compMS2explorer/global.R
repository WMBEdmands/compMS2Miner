load(file='compMS2object.RData')

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
mass.v <- sapply(object@metaData,function(x) as.numeric(x[grep("MS1_mz", names(x))][[1]][1]))
RT.v <- sapply(object@metaData, function(x) as.numeric(x[grep("MS1_RT", names(x))][[1]][1]))
TotalFeatures<-length(unique(gsub(".+_","",Features.v)))
TotalCompSpectra<-length(Features.v)