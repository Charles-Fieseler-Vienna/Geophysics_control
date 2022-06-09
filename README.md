
## Analysis of Geophysical Acoustic Emissions

## Getting started

There are two steps required:
* Run the setup_control_analysis.m file, or manually add the folders to your MATLAB path.
* Go to /src/ControlProject.m and change the "dat_foldername" property to point to the location of the data files.
  * In addition, if the subfolders are different than they appear in the dropbox, then the "dat_subfolders" array will need to be updated

## Overview of code structure

### src/
This folder includes the source code, which is split up into four folders:

* src/data
  * Source code for analysis of the data (paper_control_analysis.m) and for structuring a project (ControlProject.m)
  * Classical (non-control) analysis is in the file: paper_fft_analysis.m
  * These are organized to be run in sections marked with "%%", and will give an error if run directly
  * These will plot some examples, but NOT the paper figures
  * The main output of this section is intermediate products; See src/plot
* src/explore
  * Additional files similar to src/data, but not included in the paper 
* src/plot
  * Scripts for producing the paper figures (paper_plots.m), plus helper functions
  * This is organized to be run in sections marked with "%%", and will give an error if run directly
  * Note that additional formatting was applied after output by MATLAB to achieve the final paper version
* src/util
  * Additional utility functions

### intermediate/

Contains intermediate products produced by src/data scripts and needed by paper_plots.m

## Reproducing paper figures

All figures are generated by the script paper_plots.m

This script requires intermediate products to be produced by the paper_control_analysis.m and paper_fft_analysis.m scripts.
These intermediate products are included, but can be regenerated from raw data as obtained from the DOI listed in the paper.
