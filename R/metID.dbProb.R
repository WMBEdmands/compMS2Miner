#' selects best annotations based on substructure annotations identified
#' by \code{\link{subStructure.Annotate}} 
#' @description Most probable database annotations either automatically 
#' decided based on substructure type detected by the 
#' \code{\link{subStructure.Annotate}} or user supplied most probable annotations 
#' one composite spectrum at a time.
#' 
#' @param object a compMS2 class object
#' @param nameFeat a unique name for a composite spectrum of interest. If not 
#' supplied (default) all most probable annotations are decided automatically
#' and for all composite spectra. Previous most probable annotations will not
#' be overwritten if the function is run more than once.
#' @param DBids unique database identifier for a specific composite spectrum, in
#' combination with nameFeat argument.
#' 
#' @return a compMS2 class object with most probable annotation(s) 
#' @export
setGeneric("metID.dbProb", function(object, ...) standardGeneric("metID.dbProb"))

setMethod("metID.dbProb", signature = "CompMS2", function(object,
                                                             nameFeat = NULL,
                                                             DBids = NULL){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if(length(DBanno(object)) == 0){
    stop("the dbAnnotate function has not yet been run the CompMS2 class object does not
         contain any potential annotations")
  } else if(!is.null(nameFeat) | !is.null(DBids)){
    if(is.null(DBids)){
      stop("A character refering to a composite spectrum to must be supplied")
    }
    if(is.null(DBids)){
      stop("A character vector containing best candidates for the feature selected 
           must be supplied")
    }
    feat.indx.tmp <- which(names(BestAnno(object)) == nameFeat)
    if(length(feat.indx.tmp) == 0){
      stop("no match was found between ", nameFeat, 
           " and any composite spectrum, please check the name and try again")
    }
    if(length(BestAnno(object)) == 0){
      BestAnno.tmp <- vector("list", length(compSpectra(object)))
      BestAnno(object) <- BestAnno.tmp
      names(BestAnno(object)) <- names(compSpectra(object))
    }
    dbAnno.tmp <- DBanno(object)[[feat.indx.tmp]]
    # subset dbannotations and add in to BestAnno object slot
    dbAnno.tmp <- dbAnno.tmp[which(dbAnno.tmp$DBid %in% DBids), , drop = F]
    if(is.null(BestAnno(object)[[feat.indx.tmp]])){
      BestAnno(object)[[feat.indx.tmp]] <- dbAnno.tmp
    } else {
      dbanno.tmp <- rbind(BestAnno(object)[[feat.indx.tmp]], dbAnno.tmp)
      BestAnno(object)[[feat.indx.tmp]] <- dbAnno.tmp
    }
    
    } else {
      # automatically select the most likely candidates based on the substructures
      # identified
      if(length(BestAnno(object)) == 0){
        BestAnno.tmp <- vector("list", length(compSpectra(object)))
        BestAnno(object) <- BestAnno.tmp
        names(BestAnno(object)) <- names(compSpectra(object))
      }
      
      Best.Anno.tmp <- lapply(c(1:length(BestAnno(object))), function(x){
        SubStr.types <- unlist(compSpectra(object)[[x]][, c("Frag.ID.type", 
                                                            "interfrag.loss.type",
                                                            "Neutral.loss.type")])
        SubStr.types <- unique(c("", unlist(strsplit(SubStr.types, ";"))))
        dbAnno.tmp <- DBanno(object)[[x]]
        dbAnno.tmp <- dbAnno.tmp[dbAnno.tmp$SubStr_type %in% SubStr.types, , drop = F]
#         if(is.null(BestAnno(object)[[x]])){
          return(dbAnno.tmp)
#         } else {
#           dbAnno.tmp <- rbind(BestAnno(object)[[x]], dbAnno.tmp)
#           dbAnno.tmp <- dbAnno.tmp[duplicated(dbAnno.tmp$DBid) == F, ]
#           return(dbAnno.tmp)
#         }
      })
      names.tmp <- names(compSpectra(object))
      BestAnno(object) <- Best.Anno.tmp
      names(BestAnno(object)) <- names.tmp
      return(object)  
    }})
