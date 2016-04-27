#' CFM wrapper

setGeneric("metID.CFM", function(object, ...) standardGeneric("metID.CFM"))

setMethod("metID.CFM", signature = "CompMS2", function(object, CFMidExe = NULL, 
                                                      featureSubSet = NULL,  
                                                      keepTempFiles = F){
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("BestAnnotations function has not yet been run or no best annotations have been selected")
  } else if (all(sapply(BestAnno(object), is.null))){
    stop("BestAnnotations function has not yet been run or no best annotations have been selected")
  } else {
    if (is.null(CFMidExe)){
      tcltk::tkmessageBox(message = 'See "http://sourceforge.net/p/cfm-id/wiki/Home/" for instructions on how to install Competitive Fragmentation Modelling')
      
      CFMidExe <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{Executable} {.exe}} {{All files} *}",
                                                       title="select your CFM-id.exe file"))
      paramsLog <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{log file} {.log}} {{All files} *}",
                                                        title="select the params_log file"))
      
    }
    
    # if CFM file empty then create list
    if(length(CFModelling(object)) == 0){
      CFM.res.tmp <- vector("list", length(compSpectra(object)))
      names(CFM.res.tmp) <- names(compSpectra(object))
      CFModelling(object) <- CFM.res.tmp
    }
    
    if(is.null(featureSubSet)){
      featureSubSet <- names(compSpectra(object))[sapply(BestAnno(object), is.null) == F]
    }
    featureSubSet.indx <- featureSubSet %in% names(compSpectra(object)) 
    if(any(featureSubSet.indx == F)){
      stop("The following composite spectra subset names do not match :", 
           paste0(sapply(featsureSubSet[featureSubSet.indx == F], message), "\n"))
    }
    
    
    #     pb <- txtProgressBar(min=0,max=length(featureSubSet),style=3)
    featureSubSet <- which(names(compSpectra(object)) %in% featureSubSet)
    for(i in 1:length(featureSubSet)){
      ##progress bar
      message(length(featureSubSet) - i, " features remaining")
      flush.console()
      
      # write tmp files 
      indx.tmp <- featureSubSet[[i]]
      QueryFile <- compSpectra(object)[[indx.tmp]]
      FeatName.tmp <- names(compSpectra(object))[indx.tmp]    
      # create CFM Query spectrum
      CFMquerySpectrum <- paste0(paste(paste(QueryFile$mass, QueryFile$intensity), 
                                       collapse="\n"))
      # if keep temp files True then create feature sub-directory
      if(keepTempFiles == F){
        tmp.CFMDir <- paste0(getwd(),"/CFMresults")
        suppressWarnings(dir.create(tmp.CFMDir))
        CFMQuerySpectrumFile <- paste0(tmp.CFMDir, "/", FeatName.tmp, "_CFMquerySpectrum.txt")
      } else {
        tmp.CFMDir <- paste0(getwd(), "/CFMresults/")
        suppressWarnings(dir.create(tmp.CFMDir))
        tmp.CFMDir <- paste0(getwd(), "/CFMresults/", FeatName.tmp)
        suppressWarnings(dir.create(tmp.CFMDir))
        CFMQuerySpectrumFile <- paste0(tmp.CFMDir, "/", FeatName.tmp, "_CFMquerySpectrum.txt")
      }
      # write CFM query spectrum file   
      writeLines(CFMquerySpectrum, con = CFMQuerySpectrumFile)
      # write params file
      param_config_file <- paste0(dirname(CFMQuerySpectrumFile), "/param_config.txt")
      writeLines(paste0("ionization_mode ", 
                        ifelse(Parameters(object)$mode == "pos", 1, 2)), 
                 con = param_config_file)
      # create and write temporary sdf file
      # extract smiles codes from tmp best anno df
      bestAnno.tmp <- BestAnno(object)[[indx.tmp]]
      SMILES_tmp <- unlist(bestAnno.tmp[, grep("SMILES", colnames(bestAnno.tmp)), 
                                        drop = F])
      SMILES_tmp <- SMILES_tmp[SMILES_tmp != ""]
      
      if(nrow(bestAnno.tmp)==1) {
        #add DB ids to SMILES names
        names(SMILES_tmp) <- paste0(gsub("SMILES", "", names(SMILES_tmp)), bestAnno.tmp$DBid)
      } else {
        #add DB ids to SMILES names
        names(SMILES_tmp) <- paste0(gsub("SMILES.+", "", names(SMILES_tmp)), 
                                    bestAnno.tmp$DBid[as.numeric(gsub(".+SMILES|SMILES", "", names(SMILES_tmp)))])
      }
      SMILES_tmp <- SMILES_tmp[duplicated(names(SMILES_tmp))==F]
      
      if(length(SMILES_tmp)>0){
        # save smiles file
        message("writing local SMILES...")
        CFMquerySMI <- paste0(paste(paste(names(SMILES_tmp), SMILES_tmp), 
                                    collapse="\n"))
        Smi_file_name_tmp <- paste0(tmp.CFMDir, "/CandSmiles.txt")
        writeLines(CFMquerySMI, con = Smi_file_name_tmp )
        command <- paste0('"',CFMidExe,'" "',CFMQuerySpectrumFile,'" "', 
                          Smi_file_name_tmp, '" -1 10 0.01 "', 
                          '" "', paste0(dirname(param_config_file), "/", 
                                        basename(paramsLog)), '" "', 
                          param_config_file, '" DotProduct "',
                          gsub("\\.txt", "_results.txt", CFMQuerySpectrumFile), '"')
        log <- system(command)
        
      }
    } # end loop
    return(object)
  }
})
