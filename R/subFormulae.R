#' subtract atomic formula y from atomic formula x
#' 
#' @param x character vector of atomic formulae (must be same length as y).
#' @param y character vector of atomic formulae (must be same length as x).
subFormulae <- function(x=NULL, y=NULL){
  if(length(x) != length(y)){
    stop('x and y must be the same length.')
  }
  # split in to atoms and number x vector 
  xEleGrStr <- gsub('([[:upper:]])', ' \\1', x)
  xEleGrStr <- gsub('$', ' ', xEleGrStr)
  xEleGrStr <- gsub("^ ", '', xEleGrStr)
  
  xEle <- gsub('[[:punct:]]|[0-9]', '', xEleGrStr)
  xNEle <- gsub('[[:punct:]]|[A-z]', "", xEleGrStr)
  
  xEle <- strsplit(xEle, ' ')
  xNEle <- strsplit(xNEle, ' ') 
  # split in to atoms and number y vector 
  yEleGrStr <- gsub('([[:upper:]])', ' \\1', y)
  yEleGrStr <- gsub('$', ' ', yEleGrStr)
  yEleGrStr <- gsub("^ ", '', yEleGrStr)
  
  yEle <- gsub('[[:punct:]]|[0-9]', '', yEleGrStr)
  yNEle <- gsub('[[:punct:]]|[A-z]', "", yEleGrStr)
  
  yEle <- strsplit(yEle, ' ')
  yNEle <- strsplit(yNEle, ' ')
  # unique elements
  uniEle <- unique(c(unlist(xEle), unlist(yEle)))
  
  constAtoms <- vector('numeric', length(uniEle))
  names(constAtoms) <- uniEle
  remFormulae <- vector('character', length(x))
  # pb <- txtProgressBar(max=length(x), style = 3)
  for(i in 1:length(x)){
    # setTxtProgressBar(pb, i)
    # x formula
    xAts <- constAtoms
    xTmp <- xNEle[[i]]
    xTmp[xTmp == ''] <- 1
    xAts[xEle[[i]]] <- as.numeric(xTmp)
    # y formula
    yAts <- constAtoms
    yTmp <- yNEle[[i]]
    yTmp[yTmp == ''] <- 1
    yAts[yEle[[i]]] <- as.numeric(yTmp)
    # subtract x from y
    remAts <- xAts - yAts
    # remove zeros
    remAts <- remAts[remAts != 0]
    remAts[remAts == 1] <- ''
    # collapse to make new formula
    remFormulae[i] <- paste(paste0(names(remAts), remAts), collapse='')
  }
  return(remFormulae)
} # end function