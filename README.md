
## Analysis of Geophysical Acoustic Emissions

## Getting started

There are two steps required:
* Run the setup_Purdue_analysis.m file, or manually add the folders to your MATLAB path.
* Go to /src/PurdueProject.m and change the "dat_foldername" property to point to the location of the data files.
  * In addition, if the subfolders are different than they appear in the dropbox, then the "dat_subfolders" array will need to be updated

## Overview of code structure

### /src/
This folder includes the source code, which is split up into three types of files, not including scratch work:

* Produce figures
  * Purdue_AGU_plots.m produces all the plots that are included in the presentation
* Produce data
  * Purdue_***_analysis.m produces the intermediate data, saved in the /intermediate/ folder
* Preprocessing
  * PurdueProject.m is a class object with various properties, e.g. the data file locations


### /intermediate/
This folder includes intermediate data products, some of which take hours (on my laptop) to reproduce. The analysis functions needed to reproduce these files are not included yet.


## Running the code

### Plots
All of the plots should be reproducible via the various sections of Purdue_AGU_plots.m, once the "Getting Started" section has been completed.

### Intermediate data and analysis
The intermediate data need additional analysis functions that are not included yet. 


## Reproducing specific figures


