#' CompMS2 class
#'
#' The CompMS2 class matches MS1 data to MS2 precursors, extracts both MS1 and
#' MS2 spectral data, noise filters and annotates both MS1 features and MS2
#'  substructures.
#' 
#' @section Slots: 
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
                        network = "list",
                        DBanno = "list",
                        BestAnno = "list",
                        subStrAnno = "data.frame",
                        spectralDB = 'list',
                        Comments = "data.frame",
                        file.paths = "character",
                        Parameters = "data.frame",
                        couchDBconn = "list",
                        MetFrag = "list",
                        CFModelling = "list"))
# CMS2object <- new("CompMS2")
# CMS2object@compSpectra <- object@compSpectra
# CMS2object@metaData  <- object@metaData
# CMS2object@network  <- object@network
# CMS2object@DBanno  <- object@DBanno
# CMS2object@BestAnno  <- object@BestAnno
# CMS2object@subStrAnno  <- object@subStrAnno
# CMS2object@Comments <- object@Comments
# CMS2object@file.paths <- object@file.paths
# CMS2object@Parameters <- object@Parameters
# CMS2object@couchDBconn <- object@couchDBconn
# CMS2object@MetFrag <- object@MetFrag
# CMS2object@CFModelling <- object@CFModelling
