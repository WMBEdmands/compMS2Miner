#' composite spectra substructure annotation 
#' @param Frag_mzabs Absolute mass accuracy difference to identify neutral losses 
#' and fragments in composite spectra (default = 0.01).
#' @param SubStrs substructure data frame (default = Substructure_masses)
#' see ?Substructure_masses for details of the mandatory table fields/ format
#' @param minRelInt minimum relative intensity to consider a spectral signal
#' for substructure annotation (default = 5 i.e. 5\% rel. int.).
#' @export
setGeneric("subStructure.Annotate", function(object, ...) standardGeneric("subStructure.Annotate"))

setMethod("subStructure.Annotate", signature = "CompMS2", 
          function(object, Frag_mzabs = 0.01, SubStrs = Substructure_masses, 
                   minRelInt = 5){
  
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else {
    
    #     if(!is.null(ppm)){
    #       ppm <- Parameters(object)$Precursor.ppm
    #     }
    
    if(!all(colnames(SubStrs) %in% colnames(Substructure_masses)))
    {
      stop('column names for the substructure masses data frame supplied do not match the required
           naming structure the required column names are as follows : \n',
           paste0(1:ncol(Substructure_masses), ". ", colnames(Substructure_masses), "\n"))
    }
    # add parameters into object
    Parameters(object) <- data.frame(Parameters(object), Frag_mzabs = Frag_mzabs, 
                                     minRelInt = minRelInt)
    # mode indx
    mode.indx <- SubStrs[, Parameters(object)$mode] == 1
    # frag indx
    Fragments <- SubStrs[SubStrs$frag == 1 & mode.indx, , drop = F]
    # neut loss indx
    Neutral.losses <- SubStrs[SubStrs$Neut.loss == 1 & mode.indx, , drop = F]
    # obtain MS1 mzs from compMS2 object
    MS1_mzs <- sapply(metaData(object), function(x) unlist(x[grep("MS1_mz", names(x))])[1])
    # comp spectra
    comp_spectra.tmp <- compSpectra(object)
    #    
    #     if(Parameters(object)$nCores > 0){
    #     
    message("matching Precursor to fragment and interfragment neutral losses and fragments in ",
            length(compSpectra(object)), " composite spectra")
    flush.console()
    
    
    comp_spectra.tmp <- lapply(1:length(comp_spectra.tmp), function(x){
      # add relative intensity
      cSpectrum.tmp <- comp_spectra.tmp[[x]]
      cSpectrum.tmp$Rel_Intensity <- 100 * (cSpectrum.tmp$intensity / max(cSpectrum.tmp$intensity))
      cSpectrum.tmp <- cSpectrum.tmp[order(cSpectrum.tmp$mass), , drop = F]
      # calculate the interfragment and precursor to fragment m/z differences
      cSpectrum.tmp$interfrag.diff <- as.numeric(c(diff(cSpectrum.tmp$mass), 0))
      cSpectrum.tmp$Precursorfrag.diff <- format(as.numeric(MS1_mzs[x]  - cSpectrum.tmp$mass), scientific=F)
      Above.minRelInt <- cSpectrum.tmp$Rel_Intensity > minRelInt
      # identify fragments,  neutral losses and interfragment differences
      cSpectrum.tmp <- data.frame(cSpectrum.tmp, t(sapply(c(1:nrow(cSpectrum.tmp)), function(y){
        # indices of frag,  interfrag and NLs
        FragIndx_tmp <- which(as.numeric(cSpectrum.tmp$mass[y]) < Fragments$monoisotopic_mass+Frag_mzabs & as.numeric(cSpectrum.tmp$mass[y]) > Fragments$monoisotopic_mass-Frag_mzabs & Above.minRelInt[y])
        InterFragIndx_tmp <- which(as.numeric(cSpectrum.tmp$interfrag.diff[y]) < Neutral.losses$monoisotopic_mass+Frag_mzabs & as.numeric(cSpectrum.tmp$interfrag.diff[y]) > Neutral.losses$monoisotopic_mass-Frag_mzabs & Above.minRelInt[y])        
        NeutLossIndx_tmp <- which(as.numeric(cSpectrum.tmp$Precursorfrag.diff[y]) < Neutral.losses$monoisotopic_mass+Frag_mzabs & as.numeric(cSpectrum.tmp$Precursorfrag.diff[y]) > Neutral.losses$monoisotopic_mass-Frag_mzabs & Above.minRelInt[y])
        # collapse names and smiles of any matches
        PA_tmp <- c(Frag.ID=paste(Fragments[FragIndx_tmp, "name"], collapse=";"), 
                    Frag.ID.SMILES=paste(Fragments[FragIndx_tmp, "SMILES"], collapse=";"), 
                    interfrag.loss=paste(Neutral.losses[InterFragIndx_tmp, "name"], collapse=";"), 
                    interfrag.loss.SMILES=paste(Neutral.losses[InterFragIndx_tmp, "SMILES"], collapse=";"), 
                    Neutral.loss=paste(Neutral.losses[NeutLossIndx_tmp, "name"], collapse=";"), 
                    Neutral.loss.SMILES=paste(Neutral.losses[NeutLossIndx_tmp, "SMILES"], collapse=";"),
                    Frag.ID.type=paste(Fragments[FragIndx_tmp, "SubStructure_type"], collapse=";"), 
                    interfrag.loss.type=paste(Neutral.losses[InterFragIndx_tmp, "SubStructure_type"], collapse=";"),  
                    Neutral.loss.type=paste(Neutral.losses[NeutLossIndx_tmp, "SubStructure_type"], collapse=";"))
        
        PA_tmp <- gsub(";$", "", PA_tmp)
        return(PA_tmp)})), stringsAsFactors = F)
      return(cSpectrum.tmp)})
   
    names(comp_spectra.tmp)  <- names(compSpectra(object))
    compSpectra(object) <- comp_spectra.tmp
    return(object)
    
  } 
  })
