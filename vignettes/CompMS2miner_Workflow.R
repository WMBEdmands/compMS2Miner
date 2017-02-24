## ---- eval=FALSE---------------------------------------------------------
#  library(compMS2Miner)
#  # assign any metabolite identification comments to a new or the same "compMS2" object
#  compMS2Example_commented <- compMS2Explorer(compMS2Example)

## ---- include=FALSE------------------------------------------------------
library(knitr)
opts_knit$set(progress=FALSE)
library(compMS2Miner)

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # file path example MS1features in comma delimited csv file
#  # (see ?example_mzXML_MS1features for details).
#  MS1features_example <- system.file("extdata", "MS1features_example.csv",
#                                     package = "compMS2Miner")
#  # mzXml file examples directory (can also be .mzML or .mgf files)
#  mzXmlDir_example <- dirname(MS1features_example)
#  # observation MS1 feature table column names character vector for corrNetwork function
#  obsNames <- c(paste0(rep("ACN_80_", 6), rep(LETTERS[1:3], each=2), rep(1:2, 3)),
#                paste0(rep("MeOH_80_", 6), rep(LETTERS[1:3], each=2), rep(1:2, 3)))
#  # use parallel package to detect number of cores
#  nCores <- parallel::detectCores()
#  
#  # read in example peakTable
#  peakTable <- read.csv(MS1features_example, header=TRUE, stringsAsFactors=FALSE)
#  # create compMS2 object
#  compMS2Demo <- compMS2Construct(MS1features = peakTable,
#                                  msDataDir = mzXmlDir_example, nCores=nCores,
#                                  mode = "pos", precursorPpm = 10, ret = 20,
#                                  TICfilter = 10000)
#  
#  # View summary of compMS2 class object at any time
#  compMS2Demo

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # dynamic noise filter
#  compMS2Demo <- deconvNoise(compMS2Demo, "DNF")
#  # View summary of compMS2 class object at any time
#  compMS2Demo
#  

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # intra-spectrum ion grouping and signal summing
#  compMS2Demo <- combineMS2(compMS2Demo, "Ions")
#  # View summary of compMS2 class object at any time
#  compMS2Demo
#  #inter spectrum ion grouping and signal summing if the argument specSimFilter
#  # is supplied (values 0-1) then any spectrum below this spectral similarity
#  # threshold (dot product) will not be included in the composite spectrum generated.
#  compMS2Demo <- combineMS2(compMS2Demo, "Spectra", specSimFilter=0.8)
#  # View summary of compMS2 class object at any time
#  compMS2Demo
#  # remove any potential contaminants based on repeated isobaric precursor masses spaced by a maximum
#  # retention time gap and of high spectral similarity
#  compMS2Demo <- combineMS2(compMS2Demo, 'removeContam', maxRtGap=Inf, nContams=4)

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # annotate substructures
#  compMS2Demo <- subStructure(compMS2Demo, "Annotate")
#  # identify most probable substructure annotation based on total relative intensity
#  # explained
#  compMS2Demo <- subStructure(compMS2Demo, "prob")
#  # summary of most probable substructure annotation
#  mostProbSubStr <- subStructure(compMS2Demo, "probSummary")
#  

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # annotate any phase II metabolites detected
#  subStrMassShift <- c(42.010565, 119.004101, 176.03209, 255.988909,
#                       305.068159, 57.021464, 161.014666, 79.956817)
#  names(subStrMassShift) <- c("acetyl", "cysteine", "glucuronide",
#                              "glucuronide sulfate", "glutathione", "glycine",
#                              "mercapturate", "sulfate")
#  subStrTypesDetected <- sapply(names(subStrMassShift), function(x){
#                                grepl(x, Comments(compMS2Demo)$compound_class,
#                                      ignore.case=TRUE)})
#  subStrMassShift <- subStrMassShift[subStrTypesDetected]
#  #  common mandatory esi adducts i.e. added to those detected in the MS1 data by CAMERA
#  mandEsiAdducts <- c('[M-H2O-H]-', '[M+Na-2H]-', '[M+Cl]-', '[M+K-2H]-',	'[M-H+HCOOH]-',	
#                      "[M-H+CH3COOH+Na]-", "[M-H+CH3COOH]-")
#  # what if any substructure detected?
#  subStrMassShift <- subStrMassShift[names(subStrMassShift) %in% names(mostProbSubStr)]
#  #  common mandatory esi adducts i.e. added to those detected in the MS1 data by CAMERA
#  #  adducts supplied as structures can be automatically interpreted with the adduct2mass function
#  #  N.B. The default is to use a large number of adduct and fragments from the
#  #  paper of Stanstrup et. al.. See ?metID.dbAnnotate for details
#  mandEsiAdducts <- c('[M+H]+', '[M+NH4]+', '[M+Na]+','[M+CH3OH+H]+', '[M+CH3COONa]+',
#                      '[M+K]+')
#  # elements to consider in the database annotation. Any structure containing an
#  # element symbol not in this list will be filtered out
#  includeElements=c('C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I')
#  # annotate composite MS2 matched MS1 features to metabolomic databases (default
#  #is HMDB, also LMSD (lipidMaps), DrugBank, T3DB and ReSpect databases can also be queried).
#  #Warning: this may take 2-3 mins as large number of query masses
#  compMS2Demo <- metID(compMS2Demo, "dbAnnotate", esiAdduct=mandEsiAdducts,
#                       includeElements=includeElements)
#  

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  compMS2Demo <- metID(compMS2Demo, "dbAnnotate", metDB=LMSD, esiAdduct=mandEsiAdducts,
#                       includeElements=includeElements)

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  compMS2Demo <- metID(compMS2Demo, "dbAnnotate", metDB=drugBank, esiAdduct=mandEsiAdducts,
#                       includeElements=includeElements)

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  compMS2Demo <- metID(compMS2Demo, "dbAnnotate", metDB=T3DB, esiAdduct=mandEsiAdducts,
#                       includeElements=includeElements)

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  compMS2Demo <- metID(compMS2Demo, "dbAnnotate", metDB=ReSpect, esiAdduct=mandEsiAdducts,
#                       includeElements=includeElements)

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # match composite spectra to spectral databases as .msp files
#  compMS2Demo <- metID(compMS2Demo, 'matchSpectralDB')
#  

## ---- eval=FALSE---------------------------------------------------------
#  lbMspFiles <- paste0('https://raw.githubusercontent.com/WMBEdmands/mspFiles/master/MoNA-export-Libraries_-_LipidBlast_SMILES_', 1:5, '.msp')
#  for(i in 1:5){
#   compMS2Demo <- metID(compMS2Demo, 'matchSpectralDB', mspFile=lbMspFiles[i])
#  }
#  

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # select most probable annotations based on substructures detected
#  compMS2Demo <- metID(compMS2Demo, "dbProb", minTimesId=1)

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # predict Phase II metabolites from SMILES codes
#  compMS2Demo <- metID(compMS2Demo, "predSMILES")

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # metFrag in silico fragmentation.
#  compMS2Demo <- metID(compMS2Demo, "metFrag")

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # CFM in silico fragmentation.
#  compMS2Demo <- metID(compMS2Demo, "CFM")

## ---- collapse=TRUE, eval=FALSE------------------------------------------
#  # calculate spectral similarity network (dot product >= 0.8 default)
#  compMS2Demo <- metID(compMS2Demo, 'specSimNetwork')
#  ##############################################################
#  ### data pre-processing MS1features with MetMSLine package ###
#  ##############################################################
#  
#  if(!require(MetMSLine)){
#    # if not installed then install from github
#    devtools::install_github('WMBEdmands/MetMSLine')
#    require(MetMSLine)
#  }
#  
#  # zero fill
#  peakTable <- zeroFill(peakTable, obsNames)
#  # log transform
#  peakTable <- logTrans(peakTable, obsNames)
#  
#  ########################################################
#  ######### end MetMSLine pre-processing #################
#  ########################################################
#  
#  # add correlation network using pre-processed MS2 matched peak table
#  compMS2Demo <- metID(compMS2Demo, method='corrNetwork', peakTable, obsNames,
#                       corrMethod='pearson', corrThresh=0.90, MTC='none', MS2only=3)

## ---- eval=FALSE---------------------------------------------------------
#  # conduct mean maximum first network neighbour chemical similarity ranking
#  compMS2Demo <- metID(compMS2Demo, 'chemSim', autoPossId=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  # conduct retention prediction modelling
#  compMS2Demo <- metID(compMS2Demo, 'rtPred')
#  plot(rtPred(compMS2Demo)$rfModel)

## ---- eval=FALSE---------------------------------------------------------
#  # build consensus metabolite annotation see ?metID.buildConsensus for details
#  compMS2Demo <- metID(compMS2Demo, 'buildConsensus', autoPossId=TRUE,
#                       include=c('massAccuracy', 'inSilico', 'rtPred', 'chemSim',
#                                 'substructure'))

## ---- eval=FALSE---------------------------------------------------------
#  # use a differential evolution genetic algorithm to identify the optimum weights
#  # of each of the 6 buildConsensus parameters ()
#  compMS2Demo <- metID(compMS2Demo, 'optimConsensus', autoPossId=TRUE,
#                       include=c('massAccuracy', 'inSilico', 'rtPred', 'chemSim',
#                                 'substructure'), itermax=40)

## ---- eval=FALSE, collapse=TRUE------------------------------------------
#  # publish your app to shinyapps.io see ?publishApp for more details
#  # you may need to install the rsconnect and shinyapps packages and also sign up for a shinyapps.io account if you don't have one.
#  # quick guide here for setting up your account: http://shiny.rstudio.com/articles/shinyapps.html
#  publishApp(compMS2Demo, appName='compMS2Demo',
#             addFiles=system.file("doc", "compMS2Miner_Workflow.pdf",
#                                  package = "compMS2Miner"))
#  

## ---- eval=FALSE---------------------------------------------------------
#  publishApp(compMS2Demo, appName='compMS2Demo',
#             writeDir='C:/',
#             addFiles=system.file("doc", "compMS2Miner_Workflow.pdf",
#                                  package = "compMS2Miner"))

## ---- eval=FALSE---------------------------------------------------------
#  mzXmlFiles <- list.files(mzXmlDir_example, full.names = TRUE, pattern='\\.mzXML$')
#  publishApp(compMS2Demo, appName='compMS2Demo',
#             writeDir='C:/',
#             addFiles=c(system.file("doc", "compMS2Miner_Workflow.pdf",
#                                    package = "compMS2Miner"),
#                        system.file("extdata", "MS1features_example.csv",
#                                    package = "compMS2Miner"), mzXmlFiles))

## ---- eval=FALSE---------------------------------------------------------
#  # full path to zip file written by publishApp function
#  compMS2Demo <- compMS2Explorer('C:/compMS2Demo/compMS2Demo.zip')

## ------------------------------------------------------------------------
sessionInfo()

