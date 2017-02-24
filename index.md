[![compMS2Miner_logo](https://github.com/WMBEdmands/compMS2Miner/blob/master/inst/shiny-apps/compMS2Explorer/www/compMS2MinerLogo.png)](http://bit.ly/28QOxj6)

compMS2Miner is an R package for comprehensive and automatable annotation of metabolomic high-resolution LC-MS datasets.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.56582.svg)](http://dx.doi.org/10.5281/zenodo.56582)
latest stable release v2.2.3 (archived on zenodo).

#Purpose
A long-standing challenge of untargeted metabolomic profiling by liquid-chromatography - high resolution mass spectrometry analysis (LC-HRMS) is rapid, precise and automatable transition from unknown mass spectral features in the form of a peak-picking software output table to full metabolite identification using MS<sup>2</sup> fragmentation data.

The number of MS<sup>2</sup> spectra commonly collected in a precursor selection based experiment (often thousands in a single MS<sup>2</sup> datafile) limits the feasibility of painstaking manual interpretation of every spectrum. A degree of autonomous unknown annotation and at the very least a confident annotation of an unknowns most probable compound class is highly desirous. A holistic overview of the relationships between spectra can greatly facilitate the correct annotation of statistically relevant unknowns. When a handful of unknowns are targeted in isolation the broader context of an unknown can be easily missed and therefore putative identities poorly interpreted.

compMS2Miner is a package in the popular R programming language developed for comprehensive unknown feature annotation using peak-picker output files and MS<sup>2</sup> data files as inputs (.mzML, .mzXML, .mgf). compMS2Miner brings together many useful metabolite identification tools (see [Features](#features) section below) and is modular and every workflow method is therefore fully extensible. 

Data curation, visualization and sharing is made possible at any stage of the compMS2Miner package workflow via an application **Composite MS2 Explorer** developed with the R shiny package. The application allows the user to rapidly create their own study-specific MS<sup>2</sup> databases for each of their chromatographic methods. Additionally an msp database file can also be rapidly generated from the output of the compMS2Miner workflow.

If you find compMS2Miner useful for your metabolite annotation challenges please remember to cite us:

**compMS2Miner: an automatable metabolite identification, visualization and data-sharing R package for high-resolution LC-MS datasets**
*William Matthew Bell Edmands, Lauren M. Petrick, Dinesh Kumar Barupal, Augustin Scalbert, Mark Wilson, Jeffrey Wickliffe, and Stephen M Rappaport*
Analytical Chemistry Just Accepted Manuscript
[DOI: 10.1021/acs.analchem.6b02394](http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b02394)

#Installation

**1.** install the latest development version and all package dependencies in one-line of code directly from GitHub using the devtools package. First ensure devtools is installed, instructions can be found here: https://github.com/hadley/devtools
```{r}
devtools::install_github('WMBEdmands/compMS2Miner', dependencies=c("Depends", "Imports", "Suggests"), build_vignettes=TRUE)
```

**-Or-**

**2.** Download the latest stable release tar /zip and install package locally. The devtools install_local function ensures all package dependencies are also installed.
```{r}
# N.B. include full file path to your downloads directory
devtools::install_local('compMS2Miner_v2.2.3.tar.gz')
```

#Getting started

After compMS2Miner is installed begin by reading the package vignette *"compMS2Miner_Workflow"*
Just type ```vignette('compMS2Miner_Workflow')``` to view the pdf of the workflow. Or view an html version of the vignette on the web by clicking the image below:

[![compMS2Miner_logo](https://github.com/WMBEdmands/compMS2Miner/blob/master/inst/shiny-apps/compMS2Explorer/www/compMS2MinerLogoTutorial.png)](http://bit.ly/28T06oN)

Example data illustrating compMS2Miner is provided internal to the package and consists of a peak-picker output table from a nano-flow LC-HRMS metabolomic dataset of human blood samples and corresponding data-dependent MS<sup>2</sup> data files. An example workflow using this data is illustrated in the package vignette. 

#Features

The compMS2Miner structured workflow performs the following (v2.2.3, 2017/02/24): 
* Matches unknown mass spectral features to precursor MS<sup>2</sup> scans and constructs the "compMS2" class object. ```compMS2Construct()```
* Dynamically filters variable noise. ```deconvNoise()```: (additional noise filtration method options are in development) 
![DNF animation](https://github.com/WMBEdmands/compMS2Miner/blob/master/inst/shiny-apps/compMS2Explorer/www/DNFanimation.gif)
* Generates composite mass spectra by multiple scan signal summation and are only combined based on a minimum mean spectral similarity score. ```combineIons()```
* Possible contaminant/redundant spectra in precursor selection based MS<sup>2</sup> are removed based on the following criteria: isobaric repeats of precursor masses of a minimum retention time gap and high spectral similarity. ```combineIons.removeContam()```
* Interprets possible substructures from a literature curated database. ```subStructure()```
* Annotates unknown masses from several metabolomic databases. ```metID(method='dbAnnotate')```
* Matches spectral databases such as massbank in the NIST msp text database format. ```metID(method='matchSpectralDB')```
* Performs crude prediction of mammalian biotransformation metabolites. ```metID(method='predSMILES')```
* Calculates correlation and spectral similarity networks which can be visualized in the shiny interface. ```metID(method='corrNetwork')``` and ```metID(method='specSimNetwork')```
* Annotations for spectra connected by high correlation and/or spectral similarity are ranked based on a mean maximum chemical similarity scoring. Automatic annotation can be achieved regardless of ESI adduct/in-source fragment type. ```metID(method='chemSim')```
* Provides wrapper functions for pre-existing *in silico* fragmentation software [MetFrag](http://msbi.ipb-halle.de/MetFrag/) and [CFM](http://cfmid.wishartlab.com/). ```metID(method='metFrag')``` ```metID(method='CFM')```
* Molecular descriptor - Random Forest recursive feature elimination based Quantitative Structure Retention Relationship modelling based on the work described in [Cao et. al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419193/).
* Consensus score calculation based on 7 scoring metrics (compMS2Miner v2.2.3) namely mass accuracy, sub-structure, *in silico* fragmentation, predicted retention time deviation, mean maximum nearest neighbour chemical similarity, number of [PubMed citations](https://www.ncbi.nlm.nih.gov/pubmed) and a mean consensus score. ```metID(method='buildConsensus')```
* Differential evolution based weighted mean consensus score optimization based on minimization of the mean rank of confidently annotated spectra.```metID(method='optimConsensus')``` 
![optimConsen animation](https://github.com/WMBEdmands/compMS2Miner/blob/master/inst/shiny-apps/compMS2Explorer/www/optimConsensusAnimation.gif)
* An interactive table to record a user's decision making process or any confirmation (e.g. database entries, comments or literature DOIs for example). ```compMS2Explorer()``` or ```runGitHubApp()```
* The **Composite MS2 Explorer** app (accessible with the function ```compMS2explorer()```) can be very rapidly published to the shinyapps.io hosting site after setting up an account using the function ```publishApp()```, as a stand-alone zip file and/or as a study specific msp file which can be concatenated to generate an in-house/publically available database (```metID(method='compMS2toMsp')```).

An example **Composite MS2 Explorer** application created using the example data (within extdata of the package) is hosted on the shinyapps.io site here: 

<a href="http://bit.ly/28QOxj6" target="blank"><img src="https://github.com/WMBEdmands/compMS2Miner/blob/master/inst/shiny-apps/compMS2Explorer/www/screenshotCompMS2Explorer_260_120.png"/></a> 

Upon completion of the compMS2Miner workflow the user can then load the **Composite MS2 Explorer** app from the compMS2 object directly or from a zip file (```compMS2explorer()```) and systematically examine each composite spectrum (including any automatic annotation made) using all of the available tools provided in the  interface. Once a decision has been made on a putative annotation the user can then make potentially detailed comments in the interactive table. In this way metabolite identification decisions can be effectively and very efficiently recorded (such as links to journal articles and other pieces of evidence in support of an assignment). 

As a final step following systematic evaluation of the data presented by the **Composite MS2 Explorer** app it is intended that the user publishes the application  to the shinyapps.io site or as a self-contained zip file that can be easily viewed by others. Using the compMS2Miner function ```publishApp()``` the application can be publically deployed and explored by other investigators and all of the now read-only interactive table comments can be viewed. However, the user is still able to redeploy the application to the shinyapps.io site if any updates are necessary or just recreate the self-contained zip file. In theory, the published self-contained **Composite MS2 Explorer** app should be viewable in perpetuity.

This app publication approach could provide a feasible mechanism for transparency, generate experimental spectrum databases and a helpful way to share metabolite identification data alongside metabolomic/lipidomic publications.

Please give us your valuable feedback on anything you like/don't like and any suggestions for improvement or alternative methods you may find useful. We are always open to collaborations from fellow metabolomic investigators.

#Licence
The compMS2Miner package is licenced under the GPLv3 (http://www.gnu.org/licenses/gpl.html).
