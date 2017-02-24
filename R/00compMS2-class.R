#' compMS2-class
#'
#' @description The compMS2 class matches MS1 data to MS2 precursors, extracts 
#' both MS1 and MS2 spectral data, noise filters and annotates both MS1 features 
#' and MS2 substructures.
#' 
#' @slot compSpectra list containing the composite MS/MS spectra
#' @slot metaData list containing the combined scan header information and other
#' information for each composite spectrum.
#' @slot network list containing the correlation and spectral similarity network
#' graphs added by \code{\link{metID.corrNetwork}} and/or 
#' \code{\link{metID.specSimNetwork}}.
#' @slot DBanno list of dataframes one for each composite spectrum containing 
#' all monoisotopic mass accuracy based database matches added by 
#' \code{\link{metID.dbAnnotate}}.
#' @slot BestAnno list of dataframes one for each composite spectrum 
#' containing all monoisotopic mass accuracy based database matches taken from
#' the DBanno list subset by presence of substructure type detected by 
#' \code{\link{subStructure.Annotate}} and \code{\link{subStructure.prob}} and 
#' then added by \code{\link{metID.dbProb}}.
#' @slot subStrAnno data.frame of most likely MS/MS fragmentation substructures
#' detected added by \code{\link{subStructure.prob}}.
#' @slot spectralDB list of lists one for each composite spectrum. If a match to
#' a spectral database file has been made for a composite spectrum using 
#' \code{\link{metID.matchSpectralDB}} then within the list there are two lists,
#' \emph{dbSpectra} containing match information and \emph{entryInfo} containing
#' the database entry information for each match. 
#' @slot Comments data.frame containing all user made comments or annotations
#' automatically added using \code{\link{metID.chemSim}},
#' \code{\link{metID.buildConsensus}} or \code{\link{metID.optimConsensus}}.
#' @slot file.paths character vector with absolute path names of each MS/MS 
#' mzXML file.
#' @slot Parameters data.frame containing user supplied parameters for compMS2Miner
#' functions for reproducibility purposes.
#' @slot inSilico list containing two lists \emph{MetFrag} and \emph{CFM}. Within
#' each of these lists are contained any matches made to \emph{in silico} 
#' fragmentation data derived from the command line versions of MetFrag and CFM 
#' added by \code{\link{metID.metFrag}} and \code{\link{metID.CFM}}.
#' @slot rtPred list containing the randomForest recursive feature elimination
#' model added by \code{\link{metID.rtPred}}.
#' 
#' @section Methods:
#' \describe{
#' \item{compSpectra<-}{\code{signature(object = "compMS2")}: set \code{compSpectra} slot}
#' \item{compSpectra}{\code{signature(object = "compMS2")}: get \code{compSpectra} slot}
#' \item{BestAnno<-}{\code{signature(object = "compMS2")}: set \code{BestAnno} slot}
#' \item{BestAnno}{\code{signature(object = "compMS2")}: get \code{BestAnno} slot}
#' \item{CFM<-}{\code{signature(object = "compMS2")}: set the CFM list
#' in the \code{inSilico} slot}
#' \item{CFM}{\code{signature(object = "compMS2")}: get the CFM list
#' in the \code{inSilico} slot}
#' \item{Comments<-}{\code{signature(object = "compMS2")}: set \code{Comments} slot}
#' \item{Comments}{\code{signature(object = "compMS2")}: get \code{Comments} slot} 
#' \item{\link{combineMS2}}{\code{signature(object = "compMS2")}: various methods for
#' grouping and signal summation of fragment ions between MS/MS scans and removal of
#' possible contaminant ions.}
#' \item{\link{compMS2Explorer}}{\code{signature(object = "compMS2")}: shiny application for visualization of compMS2Miner results.}
#' \item{couchDBconn<-}{\code{signature(object = "compMS2")}: set \code{couchDBconn} slot}
#' \item{couchDBconn}{\code{signature(object = "compMS2")}: get \code{couchDBconn} slot} 
#' \item{DBanno<-}{\code{signature(object = "compMS2")}: set \code{DBanno} slot}
#' \item{DBanno}{\code{signature(object = "compMS2")}: get \code{DBanno} slot}
#' \item{\link{deconvNoise}}{\code{signature(object = "compMS2")}: ion fragment noise
#' filtration methods such as dynamic noise filtration.}
#' \item{filePaths<-}{\code{signature(object = "compMS2")}: set \code{filePaths} slot}
#' \item{filePaths}{\code{signature(object = "compMS2")}: get \code{filePaths} slot} 
#' \item{metaData<-}{\code{signature(object = "compMS2")}: set \code{metaData} slot}
#' \item{metaData}{\code{signature(object = "compMS2")}: get \code{metaData} slot}
#' \item{MetFrag<-}{\code{signature(object = "compMS2")}: set the MetFrag list
#' in the \code{inSilico} slot}
#' \item{MetFrag}{\code{signature(object = "compMS2")}: get the MetFrag list
#' in the \code{inSilico} slot}
#' \item{\link{metID}}{\code{signature(object = "compMS2")}: various methods for
#' metabolite identification including:
#' \enumerate{
#' \item monoisotopic mass based precursor MS1 to database matching with customisable
#' ESI adduct generation. \code{metID(object, method='dbAnnotate')}.
#' \item subset database annotations to \code{BestAnno(object)} slot based on
#' customizable substructures identified by \link{subStructure}. \code{metID(object, method='dbProb')}.
#' \item match to spectral database resources such as lipidBlast and massbank using
#' the NIST (National Institute of Standards and Technology) ASCII text database format based on dot product similarity and proportion of spectrum explained. \code{metID(object, method='matchSpectralDB')}.
#' \item predict Phase II biotransformation metabolites from SMILES codes 
#' contained in the \code{BestAnno(object)} slot (currently glucuronides, sulfates and glycine conjugates). \code{metID(object, method='predSMILES')}.
#' \item \emph{in silico} fragmentation prediction using the command line versions of both the MetFrag \code{metID(object, method='MetFrag')} and CFM (competitive fragmentation modelling) \code{metID(object, method='CFM')} softwares internal to compMS2Miner.
#' \item spectral similarity network calculation, both ion fragment and neutral loss
#' patterns are used to calculate dot product similarity scores \code{metID(object, method='specSimNetwork')}.
#' \item correlation network, MS1 level peak areas can be used to calculate a 
#' correlation network. \code{metID(object, method='corrNetwork')}.
#' \item 1st network (spectral similarity and/or correlation) neighbour mean maximum chemical similarity scoring (tanimoto) and automated
#' first pass metabolite identification \code{metID(object, method='chemSim')}.
#' \item retention time prediction using randomForest and molecular descriptors.\code{metID(object, method='rtPred')}.
#' \item build metabolite consensus based on optional combinations of mass accuracy,
#' spectral database matching, \emph{in silico} fragmentation matching, predicted
#' retention time, 1st network neighbour mean maximum chemical similarity and number of
#' pubMed citations. 
#' \item optimize weighted mean consensus metabolite scoring  using differential evolution algorithm. \code{metID(object, method='optimConsensus')}
#' }}
#' \item{network<-}{\code{signature(object = "compMS2")}: set \code{network} slot}
#' \item{network}{\code{signature(object = "compMS2")}: get \code{network} slot}
#'  \item{Parameters<-}{\code{signature(object = "compMS2")}: set \code{Parameters} slot}
#' \item{Parameters}{\code{signature(object = "compMS2")}: get \code{Parameters} slot} 
#' \item{publishApp}{\code{signature(object = "compMS2")}: publish compMS2 object
#'  to shinyapps.io account or as a self-contained zip file.}
#' \item{rtPred<-}{\code{signature(object = "compMS2")}: set \code{rtPred} slot}
#' \item{rtPred}{\code{signature(object = "compMS2")}: get \code{rtPred} slot}
#' \item{spectralDB<-}{\code{signature(object = "compMS2")}: set \code{spectralDB} slot}
#' \item{spectralDB}{\code{signature(object = "compMS2")}: get \code{spectralDB} slot}
#' \item{subsetCompMS2}{\code{signature(object = "compMS2")}: subset a compMS2 object
#' according to a character vector of compSpectra names. Retains all corresponding
#' slot names except networks derived from \link{corrNetwork} and \link{specSimNetwork}
#' which must be recalculated. This can be used to quickly remove any unwanted
#' spectra by name.}
#' \item{subStrAnno<-}{\code{signature(object = "compMS2")}: set \code{subStrAnno} slot}
#' \item{subStrAnno}{\code{signature(object = "compMS2")}: get \code{subStrAnno} slot} 
#' \item{\link{subStructure}}{\code{signature(object = "compMS2")}: identify
#' ion fragment substructures and neutral losses and identify the most probable
#' substructure.
#' }}
#'
#' @name compMS2-class
#' @exportClass compMS2
#' @rdname compMS2-class
#' @author WMB Edmands
#' @export
setClass("compMS2",
         representation(compSpectra = "list",
                        metaData = "list",
                        network = "list",
                        DBanno = "list",
                        # BestAnno = "list",
                        subStrAnno = "data.frame",
                        spectralDB = 'list',
                        Comments = "data.frame",
                        filePaths = "character",
                        Parameters = "data.frame",
                        couchDBconn = "list",
                        inSilico = "list", 
                        rtPred = 'list'))
# CMS2object <- new("compMS2")
# compSpectra(CMS2object) <- compSpectra(object)
# metaData(CMS2object)  <- metaData(object)
# network(CMS2object)  <- network(object)
# DBanno(CMS2object)  <- DBanno(object)
# BestAnno(CMS2object)  <- BestAnno(object)
# subStrAnno(CMS2object)  <- subStrAnno(object)
# Comments(CMS2object) <- Comments(object)
# filePaths(CMS2object) <- filePaths(object)
# Parameters(CMS2object) <- Parameters(object)
# couchDBconn(CMS2object) <- couchDBconn(object)
# MetFrag(CMS2object) <- object@inSilico$MetFrag
# CFM(CMS2object) <- object@inSilico$CFM
