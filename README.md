# beyond-tracking-climate
Data and R code for Lustenhouwer &amp; Parker, "Beyond tracking climate: niche evolution during native range expansion and its implications for novel invasions"

A preprint of the manuscript and supplementary material can be found on bioRxiv here: https://doi.org/10.1101/2021.06.09.447486

# Home directory
Contents:
- this README

- beyond-tracking-climate.Rproj
  -  R project file. Open this file in R studio to run all R scripts with the correct relative file paths.

- citations.txt
  - Reference list for all species occurrence data used in this project. 

# RScripts
Contents:
- LustenhouwerParker_start.R
  - Start here. This script will load the climate and species occurrence data needed to run subsequent scripts. 
  - Related figures and tables: Table 1

- LustenhouwerParker_PredictorSelection.R
  - This script shows how the climate variables were selected. It is stand-alone and does not need to be run for subsequent analyses.

- LustenhouwerParker_COUE.R
  - This script will run the COUE analyses (niche centroid shift, overlap, unfilling, and expansion).
  - Related figures and tables: Figure 1, Figure 3abde, Table 2 

- LustenhouwerParker_MaxEnt.R
  - This script will run the Maximum Entropy models.
  - Related figures and tables: Figure 2, Figure 3cf, Figure S1-4, Table 3

# Data
This folder contains all data files needed to reproduce the main results in the paper using the provided R scripts. Species occurrence data are available at 0.5° latitude x longitude resolution.

Contents:
- Shapefiles 
  - These consist of multiple files with the following extensions that belong together: .cpg, .dbf, .prj, .qpj, .shp, .shx
  - Australia: spatial extent of Australia, used for plotting
  - California: spatial extent of California, used for plotting
  - NativeRange: historic native range limit
  - StudyRegion: boundary of the entire study region in the Eurasian Holarctic

- Climate data rasters
  - Tif files containing the 19 WORLDCLIM variables (bio1-19) and the number of frost days for September-December (frost).
  - Data source: Climatic Research Unit (CRU) TS4.04. Harris I, Osborn TJ, Jones P, Lister D. 2020 Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. <i>Sci. Data</i> <b>7</b>, 109. (doi:10.1038/s41597-020-0453-3)
  - PastClimate.tif: Eurasia, 1901-1930
  - PresentClimate.tif: Eurasia, 1990-2019
  - PresentClimateAU.tif: Australia, 1990-2019
  - PresentClimateCA.tif: California, 1990-2019
  - PresentClimateNA.tif: North America, 1990-2019

- Species occurrence data for <i>Dittrichia graveolens</i>
  - CSV files containing a combination of occurrence and climate data.
  - Data sources: citations.txt (main folder) 
  - Variables:
    - x: decimal longitude
    - y: decimal latitude
    - bio1-bio19, frost: climate variables
    - counts.cru: placeholder column for subsequent calculations
    - species_occ: presence (1) or absence (0) of <i>Dittrichia graveolens</i>   
  - PastClimate_NativeRecords.csv: past climate data with occurrence records in the historic native range 
  - PresentClimate_AllRecords.csv: present climate data with all occurrences in the Eurasian study region
  - PresentClimate_AustraliaRecords.csv: present climate data with occurrences in Australia
  - PresentClimate_CaliforniaRecords.csv: present climate data with occurrences in California

- Target background data
  - CSV files containing occurrence data for target background species, and corresponding climate data
  - Data sources: 
    - GBIF.org. 2021 GBIF Occurrence Download for <i>Dittrichia viscosa</i> (L.) Greuter, (doi:10.15468/dl.p2mhke)
    - GBIF.org. 2021 GBIF Occurrence Download for <i>Inula</i> (L.), (doi:10.15468/dl.jj2q5f)
    - GBIF.org. 2021 GBIF Occurrence Download for <i>Pulicaria</i> Gaertn. ex Schreb., (doi: 10.15468/dl.mfjjxu) 
  - Variables: 
    - longitude (decimal)
    - latitude (decimal)
    - bio1-bio19, frost: climate variables
    - targetGBIF.counts: number of occurrences of target background species in cell
    - species_occ: presence (1) or absence (0) of target background species. Only presences are included in these files.
  - targetBgPast.csv: target background points with past climate data
  - targetBgPresent.csv: the same target background points, but with present climate data
