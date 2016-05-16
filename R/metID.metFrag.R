#' metFrag query function 
#' 
#' @description performs metFrag (msbi.ipb-halle.de/MetFrag/) insilico combinatorial
#' fragmentation. Local chemical structure data files (.sdf) are created from 
#' most probable annotation canonical SMILES codes. Temporary local sdf files 
#' and metfrag query files (.mf) are  created on a composite spectrum by 
#' composite spectrum basis and
#' insilico fragmentation performed. Results are read back into R and stored in
#' the compMS2 class object as results tables. In addition the temporary sdf,
#' mf and metfrag results files may also optionally be kept (keepTempSdf = TRUE)
#' and are saved in a subdirectory structure.
#' 
#' @param object a compMS2 class object
#' @param featureSubSet optional character vector of feature names to on which to
#' perform metFrag queries
#' otherwise the default is to perform metFrag queries on all features.
#' @param metFragJar file path to metfrag jar file the latest
#' version can be downloaded here: 
#'  \url{https://github.com/c-ruttkies/Tools/raw/master/MetFragCommandLineTool.jar}
#'  and the full path to the downloaded .jar file must be supplied. Alternatively,
#'  a file selection window (tcltk based) will appear if this argument is not supplied
#'  and the user must navigate to the jar file.
#' @param keepTempSdf logical default = FALSE, sdf, mf and results files will
#' be created temporarily otherwise temporary files will be retained in named
#' subdirectories (see details).
#' @details Results directories: for each MS1 feature matched to MS2 data a results 
#'                      directory is created. Features and their corresponding 
#'                      directories are named according to the schema ("EIC_",
#'                      MS1 feature EIC/unique ID number,"_M",mass-to-charge
#'                      ratio,"RT",retention time in seconds).
#' 
#' MetFrag query files : MetFrag query files (.mf) are saved in each result
#'                       directory. 
#'                 See http://c-ruttkies.github.io/MetFrag/projects/commandline/
#'                 for details
#'
#' Results sdf files : local sdf (chemical structure data files) are saved in 
#'                     each result directory. 
#'                     These files are used for MetFrag insilico fragmentation.
#' @return a compMS2 class object containing metFrag insilico fragmentation 
#' results.
#' @export
setGeneric("metID.metFrag", function(object, ...) standardGeneric("metID.metFrag"))

setMethod("metID.metFrag", signature = "CompMS2", function(object, 
                                                          featureSubSet=NULL,
                                                          metFragJar=NULL, 
                                                          keepTempSdf=FALSE){
  # to do list auto download metFragJar internally
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("metID.dbProb function has not yet been run or no probable annotations have been selected")
  } else if (all(sapply(BestAnno(object), is.null))){
    stop("metID.dbProb function has not yet been run or no probable annotations have been selected")
  } else if (!require(ChemmineR)){
    stop("the ChemmineR package must be installed from the Bioconductor repository to proceed.")
  } else if (!require(ChemmineOB)){
    stop("the ChemmineOB package must be installed from the Bioconductor repository to proceed.")
  } else {
    if (is.null(metFragJar)){
      tcltk::tkmessageBox(message = "Select the MetFrag .jar (Java Archive) file in the next window. This can be downloaded directly from Github at the following web address: https://github.com/c-ruttkies/Tools/raw/master/MetFragCommandLineTool.jar" )


      metFragJar <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{Java Archive file} {.jar}} {{All files} *}",
                                                         title="select your MetFrag .jar file"))
      
      # internal to package in external data
      # metFragJar <- system.file("extdata", "MetFragCommandLineTool.jar", package = "CompMS2miner")
    }
    
    # if metFrag file empty then create list
    if(length(MetFrag(object)) == 0){
      metFrag.res.tmp <- vector("list", length(compSpectra(object)))
      names(metFrag.res.tmp) <- names(compSpectra(object))
      MetFrag(object) <- metFrag.res.tmp
    }
    
    if(is.null(featureSubSet)){
      featureSubSet <- names(compSpectra(object))[sapply(BestAnno(object), is.null) == F]
    }
    featureSubSet.indx <- featureSubSet %in% names(compSpectra(object)) 
    if(any(featureSubSet.indx == F)){
      stop("The following composite spectra subset names do not match :", 
           paste0(sapply(featsureSubSet[featureSubSet.indx == F], message), "\n"))
    }
    
    
    #   pb <- txtProgressBar(min=0,max=length(featureSubSet),style=3)
    featureSubSet <- which(names(compSpectra(object)) %in% featureSubSet)
    for(i in 1:length(featureSubSet)){
      # write tmp sdf files then remove in working directory
      indx.tmp <- featureSubSet[[i]]
      QueryFile <- compSpectra(object)[[indx.tmp]]
      FeatName.tmp <- names(compSpectra(object))[indx.tmp]
      Neut.mass.tmp <- metaData(object)[[indx.tmp]]
      Neut.mass.tmp <- Neut.mass.tmp[[grep("MS1_mz", names(Neut.mass.tmp))[1]]][1]
      Neut.mass.tmp <- ifelse(Parameters(object)$mode == "pos", Neut.mass.tmp - 1.00726,
                              Neut.mass.tmp + 1.00726)
      
      # create metFrag Query File
      metFragQuery <- paste0("# Sample: ", names(compSpectra(object))[indx.tmp], "\n", 
                             "# Parent Mass: ", Neut.mass.tmp , "\n", 
                             "# Search PPM: ", Parameters(object)$precursorPpm, "\n", 
                             "# Mode:", ifelse(Parameters(object)$mode == "pos", 3, 1), "\n", 
                             "# Charge: 1\n",   
                             paste(paste(QueryFile$mass, QueryFile$intensity), collapse="\n"))
      # add MetFragQuery file location
      if(keepTempSdf == F){
        tmp.metFragDir <- paste0(getwd(),"/sdf_MetFragResults")
        suppressWarnings(dir.create(tmp.metFragDir))
        metFragQueryFile <- paste0(tmp.metFragDir, "/", FeatName.tmp, "_metFragQuery.mf")
      } else {
        tmp.metFragDir <- paste0(getwd(), "/sdf_MetFragResults/")
        suppressWarnings(dir.create(tmp.metFragDir))
        tmp.metFragDir <- paste0(getwd(), "/sdf_MetFragResults/", FeatName.tmp)
        suppressWarnings(dir.create(tmp.metFragDir))
        metFragQueryFile <- paste0(tmp.metFragDir, "/", FeatName.tmp, "_metFragQuery.mf")
      }
      # write mf query file   
      writeLines(metFragQuery, con=metFragQueryFile)
      
      # create and write temporary sdf file
      # extract smiles codes from tmp best anno df
      bestAnno.tmp <- BestAnno(object)[[indx.tmp]]
      SMILES_tmp <- unlist(bestAnno.tmp[, grep("SMILES", colnames(bestAnno.tmp)), 
                                        drop = F])
      SMILES_tmp <- SMILES_tmp[SMILES_tmp != ""]
      
      if(nrow(bestAnno.tmp) == 1) {
        #add DB ids to SMILES names
        names(SMILES_tmp) <- paste0(gsub("SMILES", "", names(SMILES_tmp)), bestAnno.tmp$DBid)
      } else {
        #add DB ids to SMILES names
        names(SMILES_tmp) <- paste0(gsub("SMILES.+", "", names(SMILES_tmp)), 
                                    bestAnno.tmp$DBid[as.numeric(gsub(".+SMILES|SMILES", "", names(SMILES_tmp)))])
      }
      SMILES_tmp <- SMILES_tmp[duplicated(names(SMILES_tmp))==F]
      
      if(length(SMILES_tmp) > 0){
        # convert first to SMILES set object and then SDF (ChemmineR/OB)
        message("converting SMILES to local SDF...")
        flush.console()
        SDF_tmp <- suppressWarnings(ChemmineR::smiles2sdf(SMILES_tmp))
        SDF_tmp <- SDF_tmp[ChemmineR::validSDF(SDF_tmp)]
        if(length(SDF_tmp) > 0){  
          # create APset atom pair descriptors files
          #       apset_tmp  <-  ChemmineR::sdf2ap(SDF_tmp)
          #assign names for identification
          # write SDF and apset to file
          
          SDF_file_name_tmp <- paste0(tmp.metFragDir, "/", FeatName.tmp, "localSDF.sdf")
          ChemmineR::write.SDF(SDF_tmp, SDF_file_name_tmp)
          
          message('sending metfrag query feature "', FeatName.tmp, '"...')
          command <- paste0('java -jar "', metFragJar,'" -D "', metFragQueryFile,'" -R "', 
                            dirname(metFragQueryFile),'" -F -d sdf -L "',
                            SDF_file_name_tmp,'"')
          log <- system(command, intern = T, show.output.on.console = F, ignore.stderr = T)
          # check if results returned
          fragFiles <- list.files(path=tmp.metFragDir, pattern='_fragments\\.sdf$')
          # if an error during system command
          if(length(fragFiles) > 0){
            sdf.results.name <- paste0(tmp.metFragDir,"/results_",names(compSpectra(object))[indx.tmp],".sdf")
            ###read back in .sdf fragments data to store in compMS2 object 
            #if(file.exists(sdf.results.name))
            #{
            sdf.results <- ChemmineR::read.SDFset(sdf.results.name)
            sdf.results.df <- as.data.frame(do.call("rbind",lapply(c(1:length(sdf.results)),function(x){
              tmp <- sdf.results[[x]]@datablock
              if(length(which(names(tmp)=="PeaksExplained"))==0)
              {
                tmp<-c(tmp,PeaksExplained="") 
              }
              return(tmp)
            })),stringsAsFactors=F)
            
            colnames(sdf.results.df)[1] <- "DBid"
            sdf.results.df$DatabaseID <- NULL 
            
            SumTotInt.tmp <- sum(compSpectra(object)[[indx.tmp]]$intensity)
            # create new mass and intensity columns in results and remove PeaksExplained
            sdf.results.df[,c("mass","intensity", "PropTotalInt")]<-t(sapply(strsplit(sdf.results.df$PeaksExplained," "),function(x){
              if(length(x)>0)
              {
                masses<-round(as.numeric(x[seq(1,length(x),2)]),digits=4)
                # add in proportion total int explained
                intensities <- round(QueryFile$intensity[round(QueryFile$mass, digits=4) %in% masses], digits = 2)
                PropTotalInt <- sum(intensities)/ SumTotInt.tmp
                intensities <- paste0(intensities, collapse="; ")
                tmp<-c(mass = paste0(masses, collapse="; "),intensity = intensities, 
                       PropTotalInt = round(PropTotalInt, digits = 2))
                return(tmp)
              } else {
                tmp <- c(mass = "",intensity = "", PropTotalInt = "")
                return(tmp)
              }
            }))
            sdf.results.df$PeaksExplained <- NULL
            sdf.results.df$SubStruc <- sdf.results.df$DBid
            sdf.results.df$DBid <- gsub(".+_", "", sdf.results.df$DBid)
            sdf.results.df$SubStruc <- ifelse(sdf.results.df$SubStruc == sdf.results.df$DBid, 
                                              "", sdf.results.df$SubStruc)
            sdf.results.df <- merge(sdf.results.df, bestAnno.tmp[, c("DBid", "DBname", "WebAddress"), drop = F])
            sdf.results.df <- sdf.results.df[order(sdf.results.df$Score, decreasing = T), , drop = F]
            sdf.results.df <- cbind(sdf.results.df$DBname, sdf.results.df$WebAddress, sdf.results.df[, -grep("DBname|WebAddress", colnames(sdf.results.df))])
            colnames(sdf.results.df)[1:2] <- c("DBname", "WebAddress")
            
            
            # add metFrag results back to object
            MetFrag(object)[[indx.tmp]] <- sdf.results.df
          
          } else {
            MetFrag(object)[[indx.tmp]] <- log
            message(log)
            flush.console()
          }
          # remove temporary files
          if(keepTempSdf == F){
            sdf.file.rem.tmp <- file.remove(list.files(tmp.metFragDir, full.names = T))
          } 
          
        }
      }
      ##progress bar
      message(length(featureSubSet) - i, " features remaining")
      flush.console()
    } # end loop
    return(object)
  }
})
