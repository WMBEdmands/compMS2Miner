## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  ## install rCharts directly from github using devtools
#  install_github("ramnathv/rCharts")

## ---- eval=FALSE---------------------------------------------------------
#  library(CompMS2miner)
#  compMS2shiny(compMS2example)

## ---- include=F----------------------------------------------------------
library(CompMS2miner)

## ---- collapse=TRUE------------------------------------------------------
# file path example MS1features in comma delimited csv file 
# (see ?example_mzXML_MS1features for details).
MS1features_example <- system.file("extdata", "MS1features_example.csv", 
                                   package = "CompMS2miner")
# mzXml file examples directory
mzXmlDir_example <- dirname(MS1features_example)
# use parallel package to detect number of cores
nSlaves <- parallel::detectCores()
# create compMS2 object
compMS2demo <- compMS2(MS1features = MS1features_example, 
                       mzXMLdir = mzXmlDir_example, nSlaves=nSlaves,
                       mode = "pos", precursorPpm = 10, ret = 10, 
                       TICfilter = 10000)
# View summary of compMS2 class object at any time 
compMS2demo

## ---- collapse=TRUE------------------------------------------------------
# dynamic noise filter
compMS2demo <- deconvNoise(compMS2demo, "DNF")
# View summary of compMS2 class object at any time 
compMS2demo 


## ---- collapse=TRUE------------------------------------------------------
# intra-spectrum ion grouping and signal summing
compMS2demo <- combineMS2(compMS2demo, "Ions")
# View summary of compMS2 class object at any time 
compMS2demo
#inter spectrum ion grouping and signal summing
compMS2demo <- combineMS2(compMS2demo, "Spectra") 
# View summary of compMS2 class object at any time 
compMS2demo


## ---- collapse=TRUE------------------------------------------------------
# annotate substructures
compMS2demo <- subStructure(compMS2demo, "Annotate")
# identify most probable substructure annotation based on total relative intensity
# explained
compMS2demo <- subStructure(compMS2demo, "prob")
# summary of most probable substructure annotation
mostProbSubStr <- subStructure(compMS2demo, "probSummary")


## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # annotate composite MS2 matched MS1 features to metabolomic databases (default
#  #is HMDB, also DrugBank, T3DB and ReSpect databases can also be queried).
#  #Warning: this may take 2-3 mins as large number of query masses
#   compMS2demo <- metID(compMS2demo, "dbAnnotate")
#  
#  # select most probable annotations based on substructures detected
#   compMS2demo <- metID(compMS2demo, "dbProb")
#  
#  # predict Phase II metabolites from SMILES codes
#   compMS2demo <- metID(compMS2demo, "predSMILES")
#  
#  # metFrag insilico fragmentation..in development
#  #compMS2demo <- metID(compMS2demo, "metFrag")
#  

