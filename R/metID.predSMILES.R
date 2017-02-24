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

setMethod("metID.predSMILES", signature = "compMS2", function(object){
  
  # error handling
  if(class(object) != "compMS2"){
    stop("argument object is not an compMS2 class object")
  } else if (length(BestAnno(object)) == 0){
    stop("No probable/ best annotations have yet been selected")
  } 
    # indx best Anno
    bestAnnoIndx <- sapply(BestAnno(object), is.null) == FALSE
   
  # 1. sulfates    
    DBanno(object)[bestAnnoIndx] <- lapply(DBanno(object)[bestAnnoIndx], function(x){
      # if none of the substructures in current Phase II list
      sulfIndx <- grepl("^sulfate$", x$SubStr_type, ignore.case = TRUE)
      if(any(sulfIndx)){
        # substr indexs
        hydroxylIndx <- lapply(gregexpr('O\\)|O$|^O', x$SMILES[sulfIndx]), function(y){
         if(y[1] != -1){
         unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
         }})
        # carboxyl
        carboxylIndx <- lapply(gregexpr('^OC\\(=O\\)|C\\(O\\)=O', x$SMILES[sulfIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # phosphoryl
        phosphoIndx <- lapply(gregexpr('OP\\(=O\\)|P\\(O\\)\\(=O\\)', x$SMILES[sulfIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # ketone 
        ketoIndx <- lapply(gregexpr('=O\\)|=O$|^O=', x$SMILES[sulfIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # amine
        amineIndx <- lapply(gregexpr('\\(N\\)', x$SMILES[sulfIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # set difference and replace functional group with conj group
        conjType <- sapply(1:sum(sulfIndx), function(y){
        remChar <- setdiff(hydroxylIndx[[y]], carboxylIndx[[y]])
        remChar <- setdiff(remChar, phosphoIndx[[y]])
        remChar <- setdiff(remChar, ketoIndx[[y]])
        if(length(remChar) == 2 & all(c(1, nchar(x$SMILES[sulfIndx][y])) %in% remChar)){
        smilesTmp <- x$SMILES[sulfIndx][y]  
        return(paste0(smilesTmp, 'S(O)(=O)=O'))  
        } else if(length(remChar) == 2){
        smilesTmp <- x$SMILES[sulfIndx][y]
        return(paste0(substr(smilesTmp, 1, remChar[1]), 'S(O)(=O)=O', substr(smilesTmp, remChar[2], nchar(smilesTmp)))) 
        } else if(length(remChar) == 1){
        smilesTmp <- x$SMILES[sulfIndx][y]  
        if(nchar(smilesTmp) == remChar){
        return(paste0(smilesTmp, 'S(O)(=O)=O'))
        } else if(substr(smilesTmp, 2, 2) == '[') {
          closeBrack <- regexpr('\\]', smilesTmp)
          return(paste0(substr(smilesTmp, 2, closeBrack), '(OS(O)(=O)=O)', substr(smilesTmp, closeBrack + 1, nchar(smilesTmp)))) 
        } else if(grepl('\\(|[0-9]', substr(smilesTmp, 2, 2))){
          return('no functional group') 
        } else {
          return(paste0(substr(smilesTmp, 2, 2), '(OS(O)(=O)=O)', substr(smilesTmp, 3, nchar(smilesTmp))))
        }
        } else if(length(remChar) > 2){
        smilesTmp <- x$SMILES[sulfIndx][y]
        firstPair <- which(diff(remChar) == 1)[1]
        remChar <- remChar[firstPair:(firstPair + 1)]
        return(paste0(substr(smilesTmp, 1, remChar[1]), 'S(O)(=O)=O', substr(smilesTmp, remChar[2], nchar(smilesTmp)))) 
        } else if(length(carboxylIndx[[y]]) > 0){
        smilesTmp <- x$SMILES[sulfIndx][y] 
        return(ifelse(carboxylIndx[[y]][1] == 1, paste0('OS(=O)(=O)', smilesTmp), paste0(substr(smilesTmp, 1, carboxylIndx[[y]][3]), 'S(O)(=O)=O', substr(smilesTmp, carboxylIndx[[y]][4], nchar(smilesTmp)))))
        } else if(length(amineIndx[[y]]) > 0){
          smilesTmp <- x$SMILES[sulfIndx][y] 
          return(paste0(substr(smilesTmp, 1, amineIndx[[y]][2]), 'S(O)(=O)=O', substr(smilesTmp, amineIndx[[y]][3], nchar(smilesTmp))))
        } else {
        return('no functional group')  
        }})
        # subset to remove any smiles with no function groups for conjugation
        x$SMILES[sulfIndx] <- conjType
        x$DBid[sulfIndx] <- paste0(x$DBid[sulfIndx], '_sulfate')
        x$SubStr_type[sulfIndx] <- ''
        # remove any rows with no funtional group detected
        x <- x[x$SMILES != 'no functional group', , drop=FALSE]
        }
        return(x)
        }) # end sulfates
  
    # 2. glucuronides    
    DBanno(object)[bestAnnoIndx] <- lapply(DBanno(object)[bestAnnoIndx], function(x){
      # if none of the substructures in current Phase II list
      glucIndx <- grepl("^glucuronide$", x$SubStr_type, ignore.case = TRUE)
      if(any(glucIndx)){
        # substr indexs
        hydroxylIndx <- lapply(gregexpr('O\\)|O$|^O', x$SMILES[glucIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # carboxyl
        carboxylIndx <- lapply(gregexpr('^OC\\(=O\\)|C\\(O\\)=O', x$SMILES[glucIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # phosphoryl
        phosphoIndx <- lapply(gregexpr('OP\\(=O\\)|P\\(O\\)\\(=O\\)', x$SMILES[glucIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # ketone 
        ketoIndx <- lapply(gregexpr('=O\\)|=O$|^O=', x$SMILES[glucIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # amine
        amineIndx <- lapply(gregexpr('\\(N\\)', x$SMILES[glucIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # set difference and replace functional group with conj group
        conjType <- sapply(1:sum(glucIndx), function(y){
          remChar <- setdiff(hydroxylIndx[[y]], carboxylIndx[[y]])
          remChar <- setdiff(remChar, phosphoIndx[[y]])
          remChar <- setdiff(remChar, ketoIndx[[y]])
          if(length(remChar) == 2 & all(c(1, nchar(x$SMILES[glucIndx][y])) %in% remChar)){
            smilesTmp <- x$SMILES[glucIndx][y]  
            return(paste0(smilesTmp, '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O'))  
          } else if(length(remChar) == 2){
            smilesTmp <- x$SMILES[glucIndx][y]
            return(paste0(substr(smilesTmp, 1, remChar[1]), '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O', substr(smilesTmp, remChar[2], nchar(smilesTmp)))) 
          } else if(length(remChar) == 1){
            smilesTmp <- x$SMILES[glucIndx][y]  
            if(nchar(smilesTmp) == remChar){
              return(paste0(smilesTmp, '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O'))
            } else if(substr(smilesTmp, 2, 2) == '[') {
              closeBrack <- regexpr('\\]', smilesTmp)
              return(paste0(substr(smilesTmp, 2, closeBrack), '(O[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O)', substr(smilesTmp, closeBrack + 1, nchar(smilesTmp)))) 
            } else if(grepl('\\(|[0-9]', substr(smilesTmp, 2, 2))){
              return('no functional group') 
            } else {
              return(paste0(substr(smilesTmp, 2, 2), '(O[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O)', substr(smilesTmp, 3, nchar(smilesTmp))))
            } 
          } else if(length(remChar) > 2){
            smilesTmp <- x$SMILES[glucIndx][y]
            firstPair <- which(diff(remChar) == 1)[1]
            remChar <- remChar[firstPair:(firstPair + 1)]
            return(paste0(substr(smilesTmp, 1, remChar[1]), '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O', substr(smilesTmp, remChar[2], nchar(smilesTmp)))) 
          } else if(length(carboxylIndx[[y]]) > 0){
            smilesTmp <- x$SMILES[glucIndx][y] 
            return(ifelse(carboxylIndx[[y]][1] == 1, paste0('O[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O', smilesTmp), paste0(substr(smilesTmp, 1, carboxylIndx[[y]][3]), '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O', substr(smilesTmp, carboxylIndx[[y]][4], nchar(smilesTmp)))))
          } else if(length(amineIndx[[y]]) > 0){
            smilesTmp <- x$SMILES[glucIndx][y] 
            return(paste0(substr(smilesTmp, 1, amineIndx[[y]][2]), '[C@@H]9O[C@@H]([C@@H](O)[C@H](O)[C@H]9O)C(O)=O', substr(smilesTmp, amineIndx[[y]][3], nchar(smilesTmp))))
          } else {
            return('no functional group')  
          }})
        # subset to remove any smiles with no function groups for conjugation
        x$SMILES[glucIndx] <- conjType
        x$DBid[glucIndx] <- paste0(x$DBid[glucIndx], '_glucuronide')
        x$SubStr_type[glucIndx] <- ''
        # remove any rows with no funtional group detected
        x <- x[x$SMILES != 'no functional group', , drop=FALSE]
      }
      return(x)
    }) # end glucuronides

    # 3. glycines
    DBanno(object)[bestAnnoIndx] <- lapply(DBanno(object)[bestAnnoIndx], function(x){
      # if none of the substructures in current Phase II list
      glyIndx <- grepl("^glycine$", x$SubStr_type, ignore.case = TRUE)
      if(any(glyIndx)){
       # carboxyl
        carboxylIndx <- lapply(gregexpr('C\\(O\\)=O', x$SMILES[glyIndx]), function(y){
          if(y[1] != -1){
            unlist(mapply(seq, from=y, to=(y + attr(y, 'match.length')) - 1, by=1, SIMPLIFY = FALSE))
          }})
        # set difference and replace functional group with conj group
        conjType <- sapply(1:sum(glyIndx), function(y){
          if(length(carboxylIndx[[y]]) > 0){
            smilesTmp <- x$SMILES[glyIndx][y] 
            return(paste0(substr(smilesTmp, 1, carboxylIndx[[y]][3]), 'NCC(O)=O', substr(smilesTmp, carboxylIndx[[y]][4], nchar(smilesTmp)))) 
          } else {
            return('no functional group')  
          }})
        # subset to remove any smiles with no function groups for conjugation
        x$SMILES[glyIndx] <- conjType
        x$DBid[glyIndx] <- paste0(x$DBid[glyIndx], '_glycine')
        x$SubStr_type[glyIndx] <- ''
        # remove any rows with no funtional group detected
        x <- x[x$SMILES != 'no functional group', , drop=FALSE]
      }
      return(x)
    }) # end glycines
    
    return(object)
}) # end function
