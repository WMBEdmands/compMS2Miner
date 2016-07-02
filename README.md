[![CompMS2miner_logo](https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/inst/shiny-apps/compMS2explorer/www/CompMS2minerLogoExApp.png)](http://bit.ly/28QOxj6)

CompMS2miner is an R package for total metabolome identification of metabolomic high-resolution LC-MS datasets.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.56582.svg)](http://dx.doi.org/10.5281/zenodo.56582)
latest stable release v1.2.5 (archived on zenodo).

#Purpose
A long-standing challenge of untargeted metabolomic profiling by liquid-chromatography - high resolution mass spectrometry analysis (LC-hrMS) is rapid, precise and automatable transition from unknown mass spectral features in the form of a peak-picking software output table to full metabolite identification.

CompMS2miner is a package in the R programming language developed for comprehensive unknown feature identification using peak-picker output files and MS/MS data files as inputs. CompMS2miner brings together many useful metabolite identification tools (see [Features](#features) section below) and is modular and therefore fully extensible. 

Data curation, visualization and sharing is made possible at any stage of the CompMS2miner package workflow via an application **Composite MS2 Explorer** developed with the R shiny package.

#Installation

**1.** install the latest development version and all package dependencies in one-line of code directly from GitHub using the devtools package. First ensure devtools is installed, instructions can be found here: https://github.com/hadley/devtools
```{r}
devtools::install_github('WMBEdmands/CompMS2miner', dependencies=c("Depends", "Imports", "Suggests"), build_vignettes=TRUE)
```

**-Or-**

**2.** You may need to install all package dependencies from CRAN and BioConductor and then download the latest stable release tar ball/zip and install package locally.
```{r}
install.packages(c('foreach', 'Rcpp', 'shiny', 'couchDB', 'fastcluster', 'data.table', 'doSNOW', 'DT', 'RcppEigen', 'reshape2', 'rjson', 'tcltk2', 'igraph', 'rhandsontable', 'rsconnect', 'shinyapps'))
 
source("https://bioconductor.org/biocLite.R")
biocLite(c('mzR', 'ChemmineR', 'ChemmineOB'))
# N.B. include full file path to your donwloads directory
install.package('CompMS2miner_v1.2.3.tar.gz', repos=NULL, type='source')
```

#Getting started

After CompMS2miner is installed begin by reading the package vignette *"CompMS2miner_Workflow"*
Just type ```vignette('CompMS2miner_Workflow')``` to view the pdf of the workflow. Or view an html version of the vignette on the web by clicking the image below:

[![CompMS2miner_logo](https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/inst/shiny-apps/compMS2explorer/www/CompMS2minerLogoTutorial.png)](http://bit.ly/28T06oN)

Example data illustrating CompMS2miner is provided internal to the package and consists of a peak-picker output table from a nano-flow LC-hrMS metabolomic dataset of human blood samples and corresponding data-dependent MS/MS data files. An example workflow using this data is illustrated in the package vignette. The CompMS2miner package is designed to offer a more complete solution to the LC-hrMS metabolite identification challenge than currently available softwares in the R language and is also complementary to other extant R packages/ workflows.

#Features

The CompMS2miner structured workflow performs the following (v1.2.5, 2016/06/28): 
* Matches unknown mass spectral features to precursor MS/MS scans and constructs the "CompMS2" class object. ```compMS2()```
* Dynamically filters variable noise. ```deconvNoise()```:
![DNF animation](https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/inst/shiny-apps/compMS2explorer/www/DNFanimation.gif)
* Generates composite mass spectra by multiple scan signal summation. ```combineIons()```
* Interprets possible substructures from a literature curated database. ```subStructure()```
* Annotates unknown masses from several metabolomic databases. ```metID(method='dbAnnotate')```
* Matches spectral databases such as massbank in the NIST msp text database format. ```metID(method='matchSpectralDB')```
* Performs crude prediction of mammalian biotransformation metabolites. ```metID(method='predSMILES')```
* Calculates correlation and spectral similarity networks which can be visualized in the shiny interface. ```metID(method='corrNetwork')``` and ```metID(method='specSimNetwork')```
* Provides wrapper functions for pre-existing *in silico* fragmentation software (http://msbi.ipb-halle.de/MetFrag/). ```metID(method='metFrag')```
* An interactive table to record a user's decision making process or any confirmation (e.g. database entries or literature DOIs). ```compMS2explorer()``` or ```runGitHubApp()```
* The **Composite MS2 Explorer** app (accessible with the function ```compMS2explorer()```) can be very rapidly published to the shinyapps.io hosting site after setting up an account using the function ```publishApp()``` or as a stand-alone zip file.

An example **Composite MS2 Explorer** application created using the example data (within extdata of the package) is hosted on the shinyapps.io site here: 

<a href="http://bit.ly/28QOxj6" target="blank"><img src="https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/inst/shiny-apps/compMS2explorer/www/screenshotCompMS2explorer_260_120.png"/></a> 

Upon completion of the CompMS2miner workflow the user can then load the **Composite MS2 Explorer** app (```compMS2explorer()```) and systematically examine each composite spectrum using all of the available tools provided in the  interface. Once a decision has been made on a putative annotation the user can then make potentially detailed comments in the interactive table. In this way metabolite identification decisions can be effectively and very efficiently recorded (such as links to journal articles and other pieces of evidence in support of an assignment). 

As a final step following systematic evaluation of the data presented by the **Composite MS2 Explorer** app it is intended that the user publishes the application  to the shinyapps.io site or as a self-contained zip file that can be easily viewed by others. Using the CompMS2miner function ```publishApp()``` the application can be publically deployed and explored by other investigators and all of the now read-only interactive table comments can be viewed. However, the user is still able to redeploy the application to the shinyapps.io site if any updates are necessary or just recreate the self-contained zip file. In theory, the published self-contained **Composite MS2 Explorer** app should be viewable in perpetuity.

This app publication approach could provide a feasible mechanism for transparency and a helpful way to share metabolite identification data alongside metabolomic/lipidomic publications.

#Licence
The CompMS2miner package is licenced under the GPLv3 (http://www.gnu.org/licenses/gpl.html).

