# Scripts for 'The causes of sea-level rise since 1900'
(c) 2020 All Rights Reserved

Authors: Thomas Frederikse, Felix Landerer, Lambert Caron, Surendra Adhikari, David Parkes, Vincent W. Humphrey, Soenke Dangendorf, Peter Hogarth, Laure Zanna, Lijing Cheng, Yun-Hao Wu

This repository contains the scripts that have been used to compute the results from the paper 'The causes of sea-level rise since 1900', Nature, 2020 (https://doi.org/10.1038/s41586-020-2591-3). The global sea-level curve and the components in NetCDF format can be obtained from NASA PO.DAAC:  https://doi.org/10.5067/GMSLT-FJPL1. A data set with the station list, VLM rates, resulting observed sea-level changes, the contributing processes, and more is available from Zenodo: https://doi.org/10.5281/zenodo.3862995. 
 

Note that this set of scripts relies on many external libraries, such as shtools (https://shtools.oca.eu/shtools/public/) for spherical harmonics, and MIDAS (http://geodesy.unr.edu/MIDAS_release/) for GPS trends.

The repository contains the following directories:

* compute_budget_terms

   The routines to compute global- and basin-mean estimates of all contributing processes. This file reads the individual gridded steric and GRD maps and averages them for each basin and the global oceans. 
* compute_grd

   All routines to prepare and generate the GRD ensemble members. This routine produces a gridded and barystatic estimate of the mass changes associated with each ensemble member. These results are read by the scripts in 'compute_budget_terms'.
* figures

   These scripts save all data in GMT-readable formats for creating the figures.
* region_selection

   Scripts used to quality-check all the individual tide-gauge and VLM observations and write the final list of tide-gauge regions used to compute basin-mean and global-mean sea-level changes. 
* results

   Routines to compute all statistics from the large ensemble. Most of the work occurs in 'compute_basin_global_stats.py'. 
* tables

   Scripts to save data in excel tables.
* tg_data

   Reads and processes all the tide-gauge data. Also contains the routines to merge the ERA5 and ERA20c reanalysis wind and pressure fields used to remove local barotropic variability from each individual tide-gauge record. 
* virtual_station

   Scripts to sample the gridded GRD fields at the tide-gauge regions and to merge the individual tide-gauge records into basin-mean and global-mean sea level using the virtual-station method.
* vlm

   All scripts to compute VLM trends from GPS and altimetry-TG records.
* GMT

   GMT scripts to plot all figures from the main text and supporting information. Get GMT here: https://github.com/GenericMappingTools/gmt. The plots use, and the scripts include some colormaps from:
  * Cynthia Brewer's ColorBrewer2 (http://colorbrewer2.org and http://soliton.vm.bytemark.co.uk/pub/cpt-city/jjg/cbcont/index.html)
  * Fabio Crameri (http://www.fabiocrameri.ch/colourmaps.php)





