[![CompMS2miner_logo](https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/CompMS2minerLogoExApp.png)](http://bit.ly/28QOxj6)

CompMS2miner is an R package for total metabolome identification of metabolomic high-resolution LC-MS datasets.

<a href="http://dx.doi.org/10.5281/zenodo.56320" target="blank"><img src="https://zenodo.org/badge/doi/10.5281/zenodo.56320.svg"/></a> 
latest stable release v1.2.3 (archived on zenodo).

Purpose
========
A long-standing challenge of untargeted metabolomic profiling by liquid-chromatography - high resolution mass spectrometry analysis (LC-hrMS) is rapid, precise and automatable transition from unknown mass spectral features in the form of a peak-picking software output peak tables to full metabolite identification.

CompMS2miner is a package in the R programming language and was developed to facilitate rapid, comprehensive unknown feature identification using peak-picker output files and MS/MS data files as inputs. 

Data curation, visualization and sharing is made possible at any stage of the CompMS2miner package workflow via an application *Composite MS2 Explorer* developed with the R shiny package.

Installation
==============

Download the latest release tar ball and install package locally

```{r}
install.package('your_downloads_directory/CompMS2miner_vLatestVersion.tar.gz', repos=NULL, type='source')
```
Or install the development version directly from github using the devtools package. First install devtools,
instructions can be found here: https://github.com/hadley/devtools

```{r}
library(devtools)

install_github('WMBEdmands/CompMS2miner', build_vignettes=TRUE)
```

Getting started
===============

After CompMS2miner is installed begin by reading the package vignette *"CompMS2miner_Workflow"*
Just type ```vignette('CompMS2miner_Workflow')``` to view the pdf of the workflow.

Example data illustrating CompMS2miner is provided consisting of a peak-picker output table of nano-flow LC-hrMS metabolomic dataset of human blood samples and data-dependent MS/MS data files, which is also made available as external example data within the CompMS2miner package. A example workflow using this data is illustrated in the package vignette. The CompMS2miner package is designed to offer a more complete solution to the LC-hrMS metabolite identification challenge than currently available softwares in the R language and is also complementary to other extant R packages/ workflows.

Features
========

The CompMS2miner structured workflow performs the following: 
* matches unknown mass spectral features to precursor MS/MS scans and constructs the "CompMS2" class object. ```compMS2()```
* dynamically filters variable noise. ```deconvNoise()```
* generates composite mass spectra by multiple scan signal summation. ```combineIons()```
* interprets possible substructures from a literature curated database. ```subStructure()```
* annotates unknown masses from several metabolomic databases. ```metID(method='dbAnnotate')```
* matches spectral databases such as massbank in the NIST msp text database format. ```metID(method='matchSpectralDB')```
* performs crude prediction of mammalian biotransformation metabolites. ```metID(method='predSMILES')```
* calculates correlation and spectral similarity networks which can be visualized in the shiny interface. ```metID(method='corrNetwork')``` and ```metID(method='specSimNetwork')```
* provides wrapper functions for pre-existing insilico fragmentation software (http://msbi.ipb-halle.de/MetFrag/). ```metID(method='metFrag')```
* an interactive table to record a user's decision making process or any confirmation (e.g. database entries or literature DOIs). ```compMS2explorer()``` or ```runGitHubApp()```
* The CompMS2explorer app can be very rapidly published to the shinyapps.io hosting site after setting up an account. ```publishApp()```

An example **Composite MS2 Explorer** application created using the example data (within extdata of the package) is hosted on the shinyapps.io site here: 

<a href="http://bit.ly/28QOxj6" target="blank"><img src="https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/screenshotCompMS2explorer_260_120.png"/></a> 

Upon completion of the CompMS2miner workflow it is intended that the user publishes the **Composite MS2 Explorer** application  to the shinyapps.io site. Using the CompMS2miner function ```publishApp``` the application can be publically deployed and explored by other invesigators and all of the now readOnly interactive table comments can be viewed. The user is still able to redeploy the application to the shinyapps.io site if necessary. 

This app publication approach could provide a feasible mechanism for transparency and a helpful way to share metabolite identification data alongside metabolomic/lipidomic publications.

licence
=============
The CompMS2miner package is licenced under the GPLv2 + the licence file (http://www.gnu.org/licenses/gpl.html).

