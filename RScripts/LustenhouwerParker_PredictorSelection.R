# R code for:
# Beyond tracking climate: niche evolution during native range expansion and its implications for novel invasions
# N. Lustenhouwer & I.M. Parker

# Script 2: Climate variable selection
# This code is not needed to run script 3 and 4; it serves to illustrate how the climate variables were selected.

# Last edit: August 27, 2021

## Run Script 1 as source code
source("Rscripts/LustenhouwerParker_start.R")

## Initial selection based on biology and redundancy ##

# temperature variability 
plot(crop(present.study.region[[c("bio4","bio2","bio7","bio3")]], extEU)) # bio4: T seasonality, bio2: mean diurnal range, bio7: T annual range, bio3: bio2/bio7
ecospat.cor.plot(present.occ[,c("bio2","bio3","bio4","bio7")])
# > keep bio2
# > between bio4 and bio7, keep bio4 

# precipitation
plot(crop(present.study.region, extEU)[["bio12"]], main="bio12") # annual precip
plot(crop(present.study.region, extEU)[["bio13"]], main="bio13") # precip of wettest month
plot(crop(present.study.region, extEU)[["bio14"]], main="bio14") # precip of driest month
plot(crop(present.study.region, extEU)[["bio15"]], main="bio15") # precip seasonality
plot(crop(present.study.region, extEU)[["bio16"]], main="bio16") # precip of wettest quarter
plot(crop(present.study.region, extEU)[["bio17"]], main="bio17") # precip of driest quarter
ecospat.cor.plot(present.occ[,c("bio12","bio13", "bio14","bio15","bio16","bio17")]) 
# > keep bio12, bio15 (correlated with bio13, bio16, bio14, bio17), and bio17

# summer and winter temperatures
plot(crop(present.study.region, extEU)[["bio5"]], main="bio5") # max T of warmest month
plot(crop(present.study.region, extEU)[["bio6"]], main="bio6") # min T of coldest month 
plot(crop(present.study.region, extEU)[["bio10"]], main="bio10") # mean T of warmest quarter
plot(crop(present.study.region, extEU)[["bio11"]], main="bio11") # mean T of coldest quarter
plot(crop(present.study.region, extEU)[["frost"]], main="frost") # number of frost days September - December
ecospat.cor.plot(present.occ[,c("bio5","bio6","bio10","bio11","frost")])
# > keep bio10 (summer T, quarter is better than month) and frost (biological relevance)

# "crossed variables": rely on both temperature and precipitation
plot(crop(present.study.region, extEU)[["bio8"]], main="bio8") # mean T of wettest Q - messy because it rains in summer in the south and in winter in the north
plot(crop(present.study.region, extEU)[["bio9"]], main="bio9") # mean T of driest Q
plot(crop(present.study.region, extEU)[["bio18"]], main="bio18") # precip of warmest Q (summer rain)
plot(crop(present.study.region, extEU)[["bio19"]], main="bio19") # precip of coldest Q (winter rain)
ecospat.cor.plot(present.occ[,c("bio8","bio9", "bio18","bio19")]) 
# > keep bio18 (summer rain), bio19 (winter rain), bio9 (how warm it is when it's dry)

# first variables selected based on biology and redundancy:
bioclimsub1 <- c("bio1","bio2","bio4","bio9","bio10","bio12","bio15","bio17","bio18","frost") 

## Final selection based on multicollinearity ##

# assess collinearity
correl <- cor(present.occ[,bioclimsub1], method="pearson") 
pos <- correl > .76
neg <- correl < -.76
spoke(pos=pos, neg=neg, lwdPos=2, lwdNeg=2, colPos="black", colNeg="red", pty="s")
vif(present.occ[,bioclimsub1]) # variance inflation factor analysis
ecospat.cor.plot(present.occ[,bioclimsub1]) # pairwise pearson correlations

# bio1 (mean annual temperature) has to go; frost can't go with bio10 or bio9
bioclimsub1b <- c("bio2","bio4","bio12","bio15","bio17","bio18","frost") 

# assess collinearity again
correl <- cor(present.occ[,bioclimsub1b], method="pearson")
pos <- correl > .76
neg <- correl < -.76
spoke(pos=pos, neg=neg, lwdPos=2, lwdNeg=2, colPos="black", colNeg="red", pty="s") # bio12 is connected to bio17 and bio18
vif(present.occ[,bioclimsub1b]) # bio12 >10
ecospat.cor.plot(present.occ[,bioclimsub1b]) # bio17 and bio12 have a pearson correlation of .81. 

# between bio12 and bio17, keep bio17 (precip in driest quarter is more biologically relevant than annual precipitation) 
bioclimsub1c <- c("bio2","bio4","bio15","bio17","bio18","frost") 

# assess collinearity again
correl <- cor(present.occ[,bioclimsub1c], method="pearson")
pos <- correl > .76
neg <- correl < -.76
spoke(pos=pos, neg=neg, lwdPos=2, lwdNeg=2, colPos="black", colNeg="red", pty="s") # looks good
vif(present.occ[,bioclimsub1c]) # max 3.9
vifstep(present.occ[,bioclimsub1c]) # no collinearity problems remaining
ecospat.cor.plot(present.occ[,bioclimsub1c]) # Matches criteria: max correlation is 0.75 (bio2~bio15)

# double-check that variables also look good for the past climate data
vif(past.occ.original[,bioclimsub1c]) # max ~4
ecospat.cor.plot(past.occ.original[,bioclimsub1c]) # max .70

# Plot final variable selection
plot(crop(past.study.region[[bioclimsub1c]], extEU))
plot(crop(present.study.region[[bioclimsub1c]], extEU))

