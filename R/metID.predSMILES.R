#' Phase II metabolite prediction
#' 
#' @description calculates possible Phase II metabolite structures from
#' canonical SMILES codes of most probable metabolite annotations. Currently
#' the algorithm predicts only certain possible Phase II metabolites and no
#' Phase I metabolism. The simple cases of acyl-, hydroxyl- and amine- sulfates
#' and glucuronides and glycine conjugates are predicted based on the presence
#' of these functional groups within the SMILES code.
#' 
#' @param object a compMS2 object (only possible when probable annotations
#' i.e. metID.dbProb has been already performed).
#' 
#' @return a compMS2 class object containing predicted Phase II metabolites from
#' most probable database annotations.
#' 
#' @export
setGeneric("metID.predSMILES", function(object, ...) standardGeneric("metID.predSMILES"))

setMethod("metID.predSMILES", signature = "CompMS2", function(object){
  
  # error handling
  if(class(object) != "CompMS2"){
    stop("argument object is not an CompMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("No probable/ best annotations have yet been selected")
  } else {
    # indx best Anno
    bestAnno.indx <- sapply(BestAnno(object), is.null) == F
    
    SubstrucStrings<-list(Hydroxyl_sulfate="S(O)(=O)=O",
                          Hydroxyl_glucuronide="[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O",
                          Acyl_sulfate="S(O)(=O)=O",
                          Acyl_glucuronide="[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O",
                          Amine_sulfate="S(O)(=O)=O",
                          Amine_glucuronide="[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O",
                          glycine="NCC(O)=O")
    
    #   bestAnno.tmp <- 
    bestAnno.tmp <- lapply(BestAnno(object)[bestAnno.indx], function(x){
      # if none of the substructures in current Phase II list
      PhaseIIsub.indx <- grep("^sulfate$|^glucuronide$|^glycine$", x$SubStr_type, ignore.case = T)
      if(length(PhaseIIsub.indx) > 0){
        # substr indx
        SMILES.sub <- x$SMILES[PhaseIIsub.indx]
        # add new pred smiles colnames
        SubstrucStringsIndx.tmp <- as.numeric(unlist(lapply(unique(x$SubStr_type[PhaseIIsub.indx]), grep, names(SubstrucStrings))))
        x[, paste0("If_", names(SubstrucStrings)[SubstrucStringsIndx.tmp], "_SMILES")] <- ""
        
        # carboxyl groups acyl-sulfates/ glucuronides and glycines
        Carboxyl.indx <- grep("C\\(O\\)=O", SMILES.sub)
        
        ###if any carboxyl groups
        if(length(Carboxyl.indx) > 0){  
          col.indx.tmp <- grep("Acyl|glycine", colnames(x))
          if(length(col.indx.tmp) > 0){
            AcylSub <- gsub("If_|_SMILES", "", colnames(x)[col.indx.tmp])
            stringSMILES <- unlist(SubstrucStrings[which(names(SubstrucStrings) %in% AcylSub)])
            x[PhaseIIsub.indx[Carboxyl.indx], col.indx.tmp] <- sapply(stringSMILES, function(y){
              predSmiles.tmp <- sub("C\\(O\\)=O", paste0("C(O", y, ")=O"), SMILES.sub[Carboxyl.indx])
            })
          }
        }
        # hydroxyl groups
        Hydroxyl.indx <- grep("\\(O\\)", SMILES.sub)
        ###if any hydroxyl groups
        if(length(Hydroxyl.indx) > 0){  
          col.indx.tmp <- grep("Hydroxyl", colnames(x))
          if(length(col.indx.tmp) > 0){
            HydroxylSub <- gsub("If_|_SMILES", "", colnames(x)[col.indx.tmp])
            stringSMILES <- unlist(SubstrucStrings[which(names(SubstrucStrings) %in% HydroxylSub)])
            x[PhaseIIsub.indx[Hydroxyl.indx], col.indx.tmp] <- sapply(stringSMILES, function(y){
              predSmiles.tmp <- sub("\\(O\\)", paste0("(O", y, ")"), SMILES.sub[Hydroxyl.indx])
              # replace potential phosphate and carboxyl
              subStringSMILES <- gsub("\\(", "\\\\(" ,y)
              subStringSMILES <- gsub("\\)", "\\\\)" , subStringSMILES)
              subStringSMILES <- gsub("\\[", "\\\\[" , subStringSMILES)
              subStringSMILES <- gsub("\\]", "\\\\]" , subStringSMILES)
              predSmiles.tmp <- gsub(paste0("P\\(O", subStringSMILES, "\\)\\(=O\\)"), "P(O)(=O)", predSmiles.tmp)
              predSmiles.tmp <- gsub(paste0("C\\(O", subStringSMILES, "\\)=O"), "C(O)=O", predSmiles.tmp)
              predSmiles.tmp <- ifelse(predSmiles.tmp == SMILES.sub[Hydroxyl.indx], "", predSmiles.tmp)
            })
          }
        }
        # amines
        Amine.indx <- grep("\\(N\\)", SMILES.sub)
        ###if any hydroxyl groups
        if(length(Amine.indx) > 0){  
          col.indx.tmp <- grep("Amine", colnames(x))
          if(length(col.indx.tmp) > 0){
            AmineSub <- gsub("If_|_SMILES", "", colnames(x)[col.indx.tmp])
            stringSMILES <- unlist(SubstrucStrings[which(names(SubstrucStrings) %in% AmineSub)])
            x[PhaseIIsub.indx[Amine.indx], col.indx.tmp] <- sapply(stringSMILES, function(y){
              predSmiles.tmp <- sub("\\(N\\)", paste0("(N", y, ")"), SMILES.sub[Amine.indx])
            })
          }
        }
        # remove any empty columns
        x <- x[, apply(x, 2, function(y) all(y == "")) == F]
      }
      return(x)
    })
    BestAnno(object)[bestAnno.indx] <- bestAnno.tmp
    return(object)
  }
})
