# Data supplement for 'The causes of sea-level rise since 1900'
(c) 2020 All Rights Reserved

Authors: Thomas Frederikse, Felix Landerer, Lambert Caron, Surendra Adhikari, David Parkes, Vincent W. Humphrey, Soenke Dangendorf, Peter Hogarth, Laure Zanna, Lijing Cheng, Yun-Hao Wu

This repository contains the scripts that have been used to compute the results from the paper. Note that this set of scripts relies on many external libraries, such as SHTNS for spherical harmonics, and MIDAS for GPS trends.

The repository contains the following directories:

* compute_budget_terms
   The routines to compute global- and basin-mean estimates of all contributing processes.
* compute_grd
   All routines to prepare and generate the GRD ensemble members.
* figures
   Scripts to save data in GMT-readable format.
* region_selection
   Scripts used to quality-check all the regions.
* results
   Routines to compute all statistics.
* tables
   Scripts to save data in excel tables.
* tg_data
   Tide-gauge data processing.
* virtual_station
   Virtual-station computations.
* vlm

   Compute VLM from GPS and altimetry-TG.
* GMT
GMT scripts to plot all figures from the main text and supporting information. Get GMT here: https://github.com/GenericMappingTools/gmt. The plots use, and the scripts include some colormaps from:
  * Cynthia Brewer's ColorBrewer2 (http://colorbrewer2.org and http://soliton.vm.bytemark.co.uk/pub/cpt-city/jjg/cbcont/index.html)
  * Fabio Crameri (http://www.fabiocrameri.ch/colourmaps.php)


