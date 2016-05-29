## ---- eval=FALSE---------------------------------------------------------
#  library(CompMS2miner)
#  compMS2explorer(compMS2example)

## ---- include=F----------------------------------------------------------
library(CompMS2miner)

## ---- collapse=TRUE------------------------------------------------------
# file path example MS1features in comma delimited csv file 
# (see ?example_mzXML_MS1features for details).
MS1features_example <- system.file("extdata", "MS1features_example.csv", 
                                   package = "CompMS2miner")
# mzXml file examples directory
mzXmlDir_example <- dirname(MS1features_example)
# observation MS1 feature table column names character vector for corrNetwork function
obsNames <- c(paste0(rep("ACN_80_", 6), rep(LETTERS[1:3], each=2), rep(1:2, 3)),
              paste0(rep("MeOH_80_", 6), rep(LETTERS[1:3], each=2), rep(1:2, 3)))
# use parallel package to detect number of cores
nCores <- parallel::detectCores()

# read in example peakTable
peakTable <- read.csv(MS1features_example, header=TRUE, stringsAsFactors=FALSE) 
# create compMS2 object
compMS2demo <- compMS2(MS1features = peakTable,
                       mzXMLdir = mzXmlDir_example, nCores=nCores,
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


## ---- collapse=TRUE------------------------------------------------------
# annotate composite MS2 matched MS1 features to metabolomic databases (default
#is HMDB, also DrugBank, T3DB and ReSpect databases can also be queried).
#Warning: this may take 2-3 mins as large number of query masses
compMS2demo <- metID(compMS2demo, "dbAnnotate")

# select most probable annotations based on substructures detected
compMS2demo <- metID(compMS2demo, "dbProb")

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # predict Phase II metabolites from SMILES codes
#  compMS2demo <- metID(compMS2demo, "predSMILES")
#  
#  # metFrag insilico fragmentation.
#  compMS2demo <- metID(compMS2demo, "metFrag")

## ---- collapse=TRUE------------------------------------------------------
##############################################################
### data pre-processing MS1features with MetMSLine package ###
##############################################################

if(!require(MetMSLine)){
  # if not installed then install from github
  devtools::install_github('WMBEdmands/MetMSLine')
  require(MetMSLine)
}

# zero fill
peakTable <- zeroFill(peakTable, obsNames)
# calculate coefficient of variation for Acetonitrile extraction replicates
peakTable <- cvCalc(peakTable, obsNames[grep('ACN', obsNames)], thresh=Inf)
# rename coeffVar column
colnames(peakTable)[ncol(peakTable)] <- 'coeffVar_ACN'
# calculate coefficient of variation for Methanol extraction replicates
peakTable <- cvCalc(peakTable, obsNames[grep('MeOH', obsNames)], thresh=Inf)
# rename coeffVar column
colnames(peakTable)[ncol(peakTable)] <- 'coeffVar_MeOH'
# all features less than 20% cv either ACN replicates or MeOH replicates
peakTable <- peakTable[peakTable$coeffVar_ACN <= 20 | peakTable$coeffVar_MeOH <= 20, ]
# deconvolute data with RamClust modified for MetMSLine
wMeanTable <- ramClustMod(peakTable, obsNames)
peakTable <- wMeanTable$wMeanPspec
# log transform
peakTable <- logTrans(peakTable, obsNames)

########################################################
######### end MetMSLine pre-processing #################
########################################################

# move eic nos to first column for corr network function
peakTable <- cbind(peakTable$EICno, peakTable)

# add correlation network using pre-processed peak table
compMS2demo <- metID(compMS2demo, method='corrNetwork', peakTable, obsNames, 
                           corrMethod='pearson', corrThresh=0.9, MTC='none')


## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # publish your app to shinyapps.io see ?publishApp for more details
#  # you may need to install the rsconnect and shinyapps packages and also sign up for a shinyapps.io account if you don't have one.
#  # quick guide here for setting up your account: http://shiny.rstudio.com/articles/shinyapps.html
#  publishApp(compMS2demo, appName='compMS2demo')
#  

