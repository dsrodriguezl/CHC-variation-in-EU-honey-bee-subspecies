# Deciphering the variation in cuticular hydrocarbon profiles of six European honey bee subspecies

This repository contains all the data and code needed to reproduce the analyses of the study titled "Deciphering the variation in cuticular hydrocarbon profiles of six European honey bee subspecies".

The repository corresponds to an R project repository that can be directly open in RStudio via the CHC-variation-in-EU-honey-bee-subspecies.Rproj file.
This makes the R project portable, as the working directory is automatically set within RStudio to the folder of the repository.

This R project uses the package renv to create a reproducible environment.
To download and install all project's dependencies, you need to use `renv::restore()` the first time you open the R project in RStudio.
This will ensure you can run the code properly to reproduce the analyses of the study.
See https://rstudio.github.io/renv/ for more details on the usage of renv package.

## Data:

samples_list.csv contains the samples metadata.

composition_data.csv contains the cuticular hydrocarbon composition data.

worldclim_tables.Rdata contains the climate data per country as downloaded from the WorldClim 2.1, via the code in the "worldclim_data.R" script.

## Code:

analysis.R is the master script. 
It contains the code to replicate the entire set of analyses of the study.
It also generates the figures for the manuscript, storing them in the folder "figs".
The code in the script relies on the CHC-variation-in-EU-honey-bee-subspecies.Rproj to determine the working directory within RStudio.

The folder "scripts" contain R scripts with code that is needed for the analysis.R script to perform the analyses of the study.



