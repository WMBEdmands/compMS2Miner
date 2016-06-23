![CompMS2miner_logo](https://raw.githubusercontent.com/WMBEdmands/CompMS2miner/master/inst/shiny-apps/compMS2explorer/www/CompMS2minerLogo.png?token=AEhLqFJ1XT2XEii5PHex0YBQQiGJ4gLDks5XdHu0wA%3D%3D)
CompMS2miner is an R package for total metabolome identification of metabolomic high-resolution LC-MS datasets.

[![DOI](https://zenodo.org/badge/21719/WMBEdmands/CompMS2miner.svg)](https://zenodo.org/badge/latestdoi/21719/WMBEdmands/CompMS2miner)

A long-standing challenge of untargeted metabolomic profiling by liquid-chromatography - high resolution mass spectrometry analysis (LC-hrMS) is rapid, precise and automatable transition from unknown mass spectral features in the form of a peak-picking software output peak tables to full metabolite identification.

CompMS2miner is a package in the R programming language and was developed to facilitate rapid, comprehensive unknown feature identification using peak-picker output files and MS/MS data files as inputs. CompMS2miner matches unknown mass spectral features to precursor MS/MS scans, dynamically filters variable noise, generates composite mass spectra by multiple scan signal summation, interprets possible substructures from a literature curated database, annotates unknown masses from several metabolomic databases, performs crude prediction of mammalian biotransformation metabolites and provides wrapper functions for pre-existing insilico fragmentation software (http://msbi.ipb-halle.de/MetFrag/).

Data curation, visualization and sharing is made possible at any stage of the CompMS2miner package workflow via an application developed with the R shiny package.

Example data illustrating CompMS2miner is provided consisting of a peak-picker output table of nano-flow LC-hrMS metabolomic dataset of human blood samples and data-dependent MS/MS data files, which is also made available as external example data within the CompMS2miner package. A example workflow using this data is available in the package vignette. The CompMS2miner package is designed to offer a more complete solution to the LC-hrMS metabolite identification challenge than currently available softwares in the R language and is also complementary to other extant R packages/ workflows.
