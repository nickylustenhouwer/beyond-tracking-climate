# R code for:
# Beyond tracking climate: niche evolution during native range expansion and its implications for novel invasions
# N. Lustenhouwer & I.M. Parker
# Submitted to Proceedings B

# Script 1: Load climate and species occurrence data

# Last edit: August 6, 2021

# Code adapted from:
  # Di Cola, V. et al. 2017. ecospat: an R package to support spatial analyses and modeling of species niches and distributions. - Ecography 40: 774–787
  # Guisan, A. 2017. Habitat suitability and distribution models, with applications in R. - Cambridge University Press
  # Smith, A. B. 2020. Best practices in species distribution modeling: a workshop in R. Available at http://www.earthskysea.org/

# Load libraries for all following scripts
library(ecospat)
library(rgeos)
library(rgdal)
library(sf)  
library(maptools)
data(wrld_simpl)
library(raster)
library(biomod2)
library(dplyr)
library(dismo)
library(maxnet)
library(rgdal)
library(rgeos)
#library(geosphere)
library(scales)
library(usdm) 

# Load packages developed by:
#Smith, A. B. 2020. Best practices in species distribution modeling: a workshop in R. Available at http://www.earthskysea.org/
#library(devtools)
#install_github('adamlilith/omnibus')
#install_github('adamlilith/statisfactory')
#install_github('adamlilith/legendary')
#install_github('adamlilith/enmSdm')
library(omnibus)
library(statisfactory)
library(legendary)
library(enmSdm)

# Color palette
# devtools::install_github("an-bui/calecopal")
library(calecopal)   


#### ========== Study region and original native range limit ============== ####
study.region <- readOGR("Data/StudyRegion.shp") # study region: Eurasian Holarctic
  plot(study.region)

native.range <- readOGR("Data/NativeRange.shp") # historic native range of Dittrichia graveolens
plot(native.range)

Australia <- readOGR("Data/Australia.shp") # outline of Australia
California <- readOGR("Data/California.shp") # outline of California
  
# spatial extent for plotting
extEU <- extent(-30,50,28,60)
  
#### =========== Climate data ============== ####

# Source: Climatic Research Unit

## Eurasia - Past ##
past.bioclim <- stack("Data/PastClimate.tif") # 1901-1930
  names(past.bioclim) <- c(paste0("bio",seq(1,19)),"frost")
past.study.region <- crop(past.bioclim, study.region) # crop to study area
past.study.region <- mask(past.study.region, study.region) 

plot(past.study.region[[1]], main = "bio1, past") 
past.climatevalues <- raster::extract(past.study.region, extent(past.study.region))   

## Eurasia - Present ##
present.bioclim <- stack("Data/PresentClimate.tif") # 1990-2019
  names(present.bioclim) <- c(paste0("bio",seq(1,19)),"frost")
present.study.region <- crop(present.bioclim, study.region) # crop to study area
present.study.region <- mask(present.study.region, study.region) 

plot(present.study.region[[1]], main = "bio1, present") 
present.climatevalues <- raster::extract(present.study.region, extent(present.study.region)) 

## North America, California, and Australia (present climate, 1990-2019) ##
present.bioclim.NA <- stack("Data/PresentClimateNA.tif") # North America 
  names(present.bioclim.NA) <- c(paste0("bio",seq(1,19)),"frost")
  plot(present.bioclim.NA[[1:6]])

present.bioclim.CA <- stack("Data/PresentClimateCA.tif") # California 
  names(present.bioclim.CA) <- c(paste0("bio",seq(1,19)),"frost")  
  plot(present.bioclim.CA[[1:6]])

present.bioclim.AU <- stack("Data/PresentClimateAU.tif") # Australia 
  names(present.bioclim.AU) <- c(paste0("bio",seq(1,19)),"frost")
  plot(present.bioclim.CA[[1:6]])

  
#### ================= Species occcurrence data =============== ####

## Eurasia - Past ##
past.occ.original <- read.csv("Data/PastClimate_NativeRecords.csv", header=T) # Past climate data with occurrence records in the historic native range
records.past <- subset(past.occ.original, species_occ==1) # subset rows with presences only
  names(records.past)[1:2] <- c("longitude", "latitude")
  nrow(records.past) # 399

plot(past.study.region[[1]], main="past climate with original native range occurrences (bio1)")
  points(latitude ~ longitude, records.past, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.2) 
  lines(native.range, col="red")  

## Eurasia - Present ##
present.occ <- read.csv("Data/PresentClimate_AllRecords.csv", header=T) # Present climate data with all occurrences
records.present <- subset(present.occ, species_occ==1) # subset rows with presences only
  names(records.present)[1:2] <- c("longitude", "latitude")
  nrow(records.present) # 746

plot(present.study.region[[1]], main="present climate with all occurrences (bio1)")
  points(latitude ~ longitude, records.present, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.2) # add expanded point
  
## California ##
calflora.occ <- read.csv("Data/PresentClimate_CaliforniaRecords.csv", header=T)
  str(calflora.occ)
calflora.records <- subset(calflora.occ, species_occ==1)
plot(present.bioclim.CA[[1]], col="lightgrey")
  points(y ~ x, calflora.records, pch=16, cex=1)

## Australia ##
Australia.occ <- read.csv("Data/PresentClimate_AustraliaRecords.csv", header=T)
  str(Australia.occ)  
Australia.records <- subset(Australia.occ, species_occ==1)
  str(Australia.records)
plot(present.bioclim.AU[[1]], col="lightgrey")
  points(y ~ x, Australia.records, pch=16, cex=.5)

which(Australia.records$y == max(Australia.records$y)) # the first record is the centroid of Australia and should be removed
Australia.records <- Australia.records[-1,]
plot(present.bioclim.AU[[1]], col="lightgrey")
  points(y ~ x, Australia.records, pch=16, cex=.5)

#### ====================== Select chosen predictors ====================== ####

predictors <- c("bio2","bio4","bio15","bio17","bio18","frost") # See script PredictorSelection for full selection procedure

# Checks for multicollinearity
ecospat.cor.plot(present.occ[,predictors]) # pairwise Pearson correlation of 0.75 or less
vif(present.occ[,predictors]) # Variance Inflation Factor < 5
  
### Table 1 ###
predictors.numeric <- c(2,4,15,17,18,20)
table1 <- data.frame("predictors" = predictors,
                     "past.min" = rep(NA,6),
                     "past.median" = rep(NA,6),
                     "past.max" = rep(NA,6),
                     "present.min" = rep(NA,6),
                     "present.median" = rep(NA,6),
                     "present.max" = rep(NA,6),
                     "NA.min" = rep(NA,6),
                     "NA.median" = rep(NA,6),
                     "NA.max" = rep(NA,6),
                     "AU.min" = rep(NA,6),
                     "AU.median" = rep(NA,6),
                     "AU.max" = rep(NA,6))

for (i in 1:6) {
  this.predictor <- predictors.numeric[i]
  
  table1[i,2] <- summary(past.bioclim[[this.predictor]])[[1]]
  table1[i,3] <- summary(past.bioclim[[this.predictor]])[[3]]
  table1[i,4] <- summary(past.bioclim[[this.predictor]])[[5]]
  
  table1[i,5] <- summary(present.bioclim[[this.predictor]])[[1]]
  table1[i,6] <- summary(present.bioclim[[this.predictor]])[[3]]
  table1[i,7] <- summary(present.bioclim[[this.predictor]])[[5]]
  
  table1[i,8] <- summary(present.bioclim.NA[[this.predictor]])[[1]]
  table1[i,9] <- summary(present.bioclim.NA[[this.predictor]])[[3]]
  table1[i,10] <- summary(present.bioclim.NA[[this.predictor]])[[5]]
  
  table1[i,11] <- summary(present.bioclim.AU[[this.predictor]])[[1]]
  table1[i,12] <- summary(present.bioclim.AU[[this.predictor]])[[3]]
  table1[i,13] <- summary(present.bioclim.AU[[this.predictor]])[[5]]
  
  print(i)
}

table1
