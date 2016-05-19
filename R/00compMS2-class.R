#' CompMS2 class
#'
#' The CompMS2 class matches MS1 data to MS2 precursors, extracts both MS1 and
#' MS2 spectral data, noise filters and annotates both MS1 features and MS2
#'  substructures.
#'  
#' This line and the next ones go into the details.
#' This line thus appears in the details as well.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{compSpectra}:}{list, containing data from compSpectra}
#'    \item{\code{slot2}:}{Object of class \code{"character"}, containing data that needs to go in slot2.}
#'  }
#' @name CompMS2 
#' @rdname CompMS2
#' @aliases CompMS2-class
#' @exportClass CompMS2
#' @author WMB Edmands
setClass("CompMS2",
         representation(compSpectra = "list",
                        metaData = "list",
                        # MS1features = "data.frame",
                        DBanno = "list",
                        BestAnno = "list",
                        subStrAnno = "data.frame",
                        Comments = "list",
                        file.paths = "character",
                        Parameters = "data.frame",
                        couchDBconn = "list",
                        MetFrag = "list",
                        CFModelling = "list"))
#CMS2object <- new("CompMS2")