#' Annotate unknown features in a CompMS2 class object to database 
#' entries based on monoisotopic mass
#' 
#' @description unknown metabolite identification. MS1 features within a CompMS2
#' object are matched against a metabolite database based on ppm error tolerance.
#' Possible electrospray adducts and substructure mass shifts are calculated for
#' all data base entries. Substructure mass shifts are supplied in a 
#' substructure masses table (see default Substructure_masses).
#' 
#' @param object a compMS2 class object.
#' @param ppm a ppm error to match MS1 features to data base entries. 
#' @param metDB a metabolite data base (see default ?HMDB for table format).
#' Other currently available data bases include DrugBank, T3DB and ReSpect 
#' databases. Matching using multiple databases is also possible. 
#' @param SubStrs substructures table containing substructure mass shifts 
#' (see default ?Substructure_masses for table format)   
#' @return a compMS2 class object containing potential metabolite annotations.
#' @source \url{http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/}
#' @export
setGeneric("metID.dbAnnotate", function(object, ...) standardGeneric("metID.dbAnnotate"))

setMethod("metID.dbAnnotate", signature = "CompMS2", function(object, ppm = NULL, 
                                                        metDB = HMDB, #pH = 2,
                                                        SubStrs = Substructure_masses
){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    if(is.null(ppm)){
      ppm <- Parameters(object)$precursorPpm
    }
    
    if(!all(c("Unique_DB_ID", "monoisotopic_weight", "name", "SMILES") %in% colnames(metDB)))
    {
      stop('metDB colnames for unique database IDs, monoisotopic weights, names and canonical SMILES of data base
           entries must be named "Unique_DB_ID", "monoisotopic_weight", "name" and "SMILES" respectively' )
    }
    # obtain MS1 mzs from compMS2 object
    unknowns <- sapply(metaData(object), function(x) unlist(x[grep("MS1_mz", names(x))])[1])
    # DB entry masses
    metMasses <- metDB$monoisotopic_weight
    # DB unique IDs
    DB_unique_IDs <- metDB$Unique_DB_ID
    # DB unique names
    DB_entry_names <- metDB$name
    DB_SMILES <- metDB$SMILES
    # web address
    dbEntryWebAddress <- metDB$WebAddress
    # tmp substr mode indx
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
    
    
    DBAnnoMatches <- monoMassMatch(unknowns = unknowns, metMasses.df =  metDB,  
                                   subStrMasses = subStrMasses, subStrNames = subStrNames, 
                                   nCores = Parameters(object)$nCores, mode = Parameters(object)$mode, 
                                   ppm = Parameters(object)$precursorPpm,
                                   dbEntryWebAddress = dbEntryWebAddress)
    if(length(DBanno(object)) > 0){
      #rbind previous DB search results
      DBAnnoMatches <- lapply(c(1:length(DBanno(object))), function(x){
        tmp.df <- rbind(DBanno(object)[[x]], DBAnnoMatches[[x]])
        # remove duplicates
        tmp.df <- tmp.df[duplicated(tmp.df$DBid) == F, ]
        return(tmp.df)
      })
    } 
    names(DBAnnoMatches) <- names(compSpectra(object))
    DBanno(object) <- DBAnnoMatches
    return(object)
    
    } 
})
