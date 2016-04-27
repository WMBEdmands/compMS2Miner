#' Create cytoscape map files based on correlation (Spearman's rho)
#' @export
setGeneric("cytoScapeMap", function(object, ...) standardGeneric("cytoScapeMap"))

setMethod("cytoScapeMap", signature = "CompMS2", function(object, sampleIDstring=NULL,
                                                          corrThresh=0.7, ...) {
  
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (is.null(sampleIDstring)){
    stop("argument sampleIDstring is missing with no default")
  } else {
    message("Creating interfeature correlation matrix...")
    flush.console()
    # eic numbers from compSpectra names
    specEics <- as.numeric(gsub(".+_| ", "", names(compSpectra(object))))
    MSfeat.df <- MS1features(object)
    # subset feat.df
    logiIndx.tmp <- MSfeat.df[, 1] %in% specEics
    feat.df <- MSfeat.df[logiIndx.tmp, grep(sampleIDstring, 
                                                        colnames(MSfeat.df)), drop=F]
    # calc correlation matrix
    message("Calculating correlation matrix for ", nrow(feat.df), " features")
    flush.console()
    cor.m <- cor(t(feat.df), method="spearman")
    colnames(cor.m)  <- MSfeat.df[logiIndx.tmp, 1]
    row.names(cor.m) <- MSfeat.df[logiIndx.tmp, 1]
    # replace upper tri with zero
    cor.m[upper.tri(cor.m, diag=T) == T] <- 0
    # ID features above below corrThresh
    sif <- apply(cor.m, 2, function(x){
                                         pos.tmp <- which(x > corrThresh)
                                         neg.tmp <- which(x < -corrThresh)
                                         return(list(pos=pos.tmp, neg=neg.tmp))
                                        })
    # melt list result
    sif.df <- reshape2::melt(sif)
     
    # create sif file names
    sif.df[, 4] <- names(sif)[sif.df[, 1]]
    
    # n nodes and edges
    nodeAttr <- data.frame(name=as.numeric(unique(c(sif.df[, 3], sif.df[, 4]))))
    message(nrow(nodeAttr), 
            " nodes with ", nrow(sif.df), " edges identified at a corrThresh of ", 
            corrThresh)
    flush.console()
    # write sif file
     res.dir <- paste0(dirname(filePaths(object)[1]), 
                      "/CompMSminer_cytoScape_", corrThresh)
    message("writing results output sif, node and edge attribute files to mzXML file directory :\n", 
    dirname(res.dir))
    flush.console()
    # write sif
    write.table(paste0(sif.df[, 3], " \t tm \t ", sif.df[, 4], "\n"), 
                paste0(res.dir, ".sif"), col.names=F, quote=F, row.names=F)
    # create node names
    sif.df[, 5] <- paste0(sif.df[, 3], " (tm) ", sif.df[, 4])
    
    edgeAttr <- data.frame(name=sif.df[, 5], correl_direction=sif.df[, 2])
    # write edge attributes
    write.table(edgeAttr,
                paste0(res.dir,"_edgeAttr.txt"), sep="\t", row.names=F, quote=F)
    # create custom node attributes table
    # extract relevant info from object information
    message("creating custom node attributes table...")
    flush.console()
    compMS2_nodeAttr <- t(sapply(1:length(specEics), function(x){
      mData.tmp <- metaData(object)[[x]]
      names(mData.tmp) <- gsub(".+_", "", names(mData.tmp))
      mData.tmp <- mData.tmp[which(duplicated(names(mData.tmp)) == F)]
      mData.tmp <- unlist(mData.tmp, recursive=T)
      names(mData.tmp) <- gsub("[0-9]$", "", names(mData.tmp))
      mData.tmp <- mData.tmp[which(duplicated(names(mData.tmp)) == F)]
      return(mData.tmp)
    }))
    
    colnames(compMS2_nodeAttr)[grep("EICno", colnames(compMS2_nodeAttr))] <- "name"
    # if substructure annotations prob carried out then add info
    subStrAnno.df <- subStrAnno(object)
    
    if(nrow(subStrAnno.df) > 0){
      subStrAnno.df$SumRelInt[which(subStrAnno.df$SumRelInt == "no substructure detected")] <- 0
      subStrAnno.df$SumRelInt <- round(as.numeric(subStrAnno.df$SumRelInt), digits=2)
      name.f <- as.factor(gsub(".+_", "", subStrAnno.df[, "compSpecName"]))
      subStrAnno.by <- by(subStrAnno.df, name.f, function(x){
                c(paste(x$SubStrType, collapse="; "),
                  paste(x$SumRelInt, collapse="; "),
                  paste(x$Freq, collapse="; "))
                })
      subStrAnno.df <- do.call("rbind", subStrAnno.by)
      subStrAnno.df <- data.frame(name=as.numeric(names(subStrAnno.by)),
                                  subStrAnno.df, stringsAsFactors=F)
      colnames(subStrAnno.df)[2:4] <- c("SubStrTypes", "SumRelInt", "Freq")
      # identify largest number of substr types
      sStr.split <- strsplit(subStrAnno.df$SubStrTypes, "; ")
      sStr.length <- sapply(sStr.split, length)
      maxSstr.type <- max(sStr.length)
      subStrAnno.df <- cbind(matrix("", ncol=maxSstr.type, 
                                    nrow=nrow(subStrAnno.df)), subStrAnno.df, 
                             stringsAsFactors=F)
      colnames(subStrAnno.df)[1:maxSstr.type] <- paste0("SubStrType_", 1:maxSstr.type)
      # add subStr types 
      for(i in 1:nrow(subStrAnno.df)){
        subStrAnno.df[i, 1:(0 + sStr.length[i])] <- sStr.split[[i]]
      }
      # identify if any evidence of a substr type
      AllsStr <- unique(unlist(sStr.split))
      AllsStr <- AllsStr[AllsStr != ""]
      # create empty columns
      subStrAnno.df[, AllsStr] <- 0
      
      message(length(AllsStr), " different substructure types were detected...")
      flush.console()
      
      for(i in 1:length(AllsStr)){
        subStrAnno.df[, AllsStr[i]] <- apply(subStrAnno.df[, 1:maxSstr.type], 1,
                                             function(x) any(x == AllsStr[i]))
      }
        # boolean any Phase II metabolites
      PhaseIImetabs <- c("diglucuronide", "cysteine", "glutathione", 
                         "mercapturate", "sulfate", "glucuronide", 
                         "glucuronide sulfate", "glycine", "acetyl")
      subStrAnno.df$PhaseII <- apply(subStrAnno.df[, PhaseIImetabs], 1,
                                        function(x) any(x == T))
      #most likely substructure column
      subStrAnno.df$MostLikely_SubStructure <- subStrAnno.df$SubStrType_1
      # bind metadata to subStr info
      compMS2_nodeAttr <- merge(subStrAnno.df, compMS2_nodeAttr, by="name")
      
    }
    # merge nodeAttr name to compMS2_nodeAttr
    nodeAttr <- merge(nodeAttr, compMS2_nodeAttr, by="name")
   write.table(nodeAttr,
                paste0(res.dir,"_nodeAttr.txt"), sep="\t", row.names=F, quote=F)
    
    message("...DONE")
    flush.console()
  }
})