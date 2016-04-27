#' example mzXML files and MS1 feature table (subset to 820 -- 940 seconds)
#' 
#' @description MS1features_example.csv 3720 MS1 features from XCMS diffreport peak table
#' from a study comparing repeat extractions of human dried blood spot samples
#' using 80\% acetonitrile (ACN) to 80\% methanol (MeOH).
#' Both extraction solvents consist of repeat preparations of the same 
#' sample (A, B, C) and repeat injections (1, 2) of each preparation 
#' (i.e. A1, A2, B1, B2, C1, C2).
#' The variables are as follows :
#' \itemize{
#'  \item EICno XCMS extracted ion chromatograms from XCMS peak tables 
#'  (62 -- 24328).
#'  \item mzmed median mass-to-charge (71.0853 -- 999.6138) 
#'  \item rtmed retention time in seconds (820.01 -- 939.907)
#'  \item ACN_80_A1 80% acetonitrile extract peak areas prep. replicate A inj. 1
#'  \item ACN_80_A2	80% acetonitrile extract peak areas prep. replicate A inj. 2
#'  \item ACN_80_B1	80% acetonitrile extract peak areas prep. replicate B inj. 1
#'  \item ACN_80_B2	80% acetonitrile extract peak areas prep. replicate B inj. 2
#'  \item ACN_80_C1	80% acetonitrile extract peak areas prep. replicate C inj. 1
#'  \item ACN_80_C2	80% acetonitrile extract peak areas prep. replicate C inj. 2
#'  \item MeOH_80_A1 80% methanol extract peak areas prep. replicate A inj. 1	
#'  \item MeOH_80_A2	80% methanol extract peak areas prep. replicate A inj. 2
#'  \item MeOH_80_B1	80% methanol extract peak areas prep. replicate B inj. 1
#'  \item MeOH_80_B2	80% methanol extract peak areas prep. replicate B inj. 2
#'  \item MeOH_80_C1	80% methanol extract peak areas prep. replicate C inj. 1
#'  \item MeOH_80_C2 80% methanol extract peak areas prep. replicate C inj. 2
#'  }
#' @docType data
#' @keywords datasets
#' @name example_mzXML_MS1features
#' @format A comma delimited text file with 3720 rows and 15 variables and two
#' data-dependent MS2 files in centroid mode converted to the mzXML open format 
#' using MSConvert software (ProteoWizard 3.0.6965 64 bit) for each extraction 
#' type.
#' Data were acquired on an Agilent 6550 q-tof interfaced with a nano-flow chip 
#' cube running a small molecule C18-chip.  
NULL