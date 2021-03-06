# R code for:
# Beyond tracking climate: niche evolution during native range expansion and its implications for novel invasions
# N. Lustenhouwer & I.M. Parker

# Script 4: MaxEnt analysis

# Last edit: April 8, 2022

# Code adapted from:
  # Di Cola, V. et al. 2017. ecospat: an R package to support spatial analyses and modeling of species niches and distributions. - Ecography 40: 774–787
  # Guisan, A. 2017. Habitat suitability and distribution models, with applications in R. - Cambridge University Press
  # Smith, A. B. 2020. Best practices in species distribution modeling: a workshop in R. Available at http://www.earthskysea.org/

#### ================ Load data ========================================

## Load climate and occurrence data from starting script 
source("RScripts/LustenhouwerParker_start.R")
head(records.past)
head(records.present) 

## Load target background points 
targetBg.raw.past <- read.csv("Data/targetBgPast.csv", header=T)
  str(targetBg.raw.past)
targetBg.raw.present <- read.csv("Data/targetBgPresent.csv", header=T)
  str(targetBg.raw.present)

## Plot all presences and bg points (Figure S1)

#png("FigureS1.png", width=5, height=5, unit="in", res=300) # run to save the figure
par(mfrow=c(2,1), mar=c(2, 3, 1, 0))
plot(present.study.region[[1]], col="lightgrey", legend=F, main="")
  #lines(native.range)
  points(latitude ~ longitude, records.present, pch=15, cex=.08, col="darkblue") # expanded native range
plot(present.study.region[[1]], col="lightgrey", legend=F, main="")
  #lines(native.range)
  points(latitude ~ longitude, targetBg.raw.past, pch=15, col="darkblue", cex=.08) # target background points
dev.off()
  
plot(crop(present.study.region[[1]], extEU), col="lightgrey", legend=F, main="Presence data, zoomed in to Europe") # zoomed in to Europe
  lines(native.range)
  points(latitude ~ longitude, records.present, pch=16, cex=.5, col="white") # expanded native range
  points(latitude ~ longitude, records.past, pch=16, cex=.5) # historic native range

## Make training data frame with predictors and vector of 1/0 for presence/background
trainData.raw.past <- rbind(records.past[ , predictors], targetBg.raw.past[ , predictors])
presBg.raw.past <- c(rep(1, nrow(records.past)), rep(0, nrow(targetBg.raw.past)))
length(presBg.raw.past)

trainData.raw.present <- rbind(records.present[ , predictors], targetBg.raw.present[ , predictors])
presBg.raw.present <- c(rep(1, nrow(records.present)), rep(0, nrow(targetBg.raw.present)))
length(presBg.raw.present)


#### ================= MaxEnt =========================================================

## Model fitting 

# Feature selection and beta/theta regularization
#trainData.raw.past.combined <- cbind(presBg.raw.past, trainData.raw.past)
#tunedModel.raw.past <- trainMaxNet(data=trainData.raw.past.combined,
#                                   regMult=c(0.5, 1, 2, 3, 4, 5, 7, 8, 9, 10), 
#                                   verbose=TRUE)

# Past Model, fast workaround without tuning (using the feature classes and regularization multiplier that come out of line 55-59)
f <- maxnet.formula(p=as.vector(presBg.raw.past), data=trainData.raw.past, classes='lpq') 
tunedModel.raw.past <- maxnet(p=presBg.raw.past, data=trainData.raw.past, f=f, regmult=0.5)

# Present Model, fit with the same settings as the tuned Past Model
f <- maxnet.formula(p=as.vector(presBg.raw.present), data=trainData.raw.present, classes='lpq')  
Model.raw.present <- maxnet(p=presBg.raw.present, data=trainData.raw.present, f=f, regmult=0.5)


## Response functions for each predictor (Figure S3 and S4) ##

# Past Model #
# get min/max value of each predictor across study region
minPred.past <- minValue(past.study.region)
maxPred.past <- maxValue(past.study.region)
names(minPred.past) <- names(maxPred.past) <- names(past.study.region)

# get median value of each predictor across species' thinned presences
medianPred.past <- apply(records.past[ , predictors], 2, median)
medianPred.past

# make data frame with median value of each predictor
env.past <- as.data.frame(medianPred.past)
env.past <- t(env.past)
env.past <- env.past[rep(1, 100), ]
row.names(env.past) <- 1:nrow(env.past)
head(env.past)

par(mfrow=c(2,3), mar=c(4,4,1,1), ps=8, family="sans")
  for (pred in predictors) {
  # make copy of data frame
  thisEnv <- env.past
  # now vary focal predictor from min to max value... 
  # all other predictors keep median value
  thisEnv[ , pred] <- seq(minPred.past[pred], maxPred.past[pred], length.out=100)
  # make prediction using this data frame
  prediction <- predict(tunedModel.raw.past, thisEnv, type='cloglog')
  # plot
  plot(x=thisEnv[ , pred], y=prediction, ylim=c(0, 1), xlab=pred,
       ylab='Suitability', main="", type='l', col='black',
       lty=1, lwd=1, cex.lab=1.3)
  # add species' presences (top rug)
  rug(records.past[ , pred], side=3, col='darkblue')
  # add background sites (bottom rug)
  rug(targetBg.raw.past[ , pred], side=1, col='darkblue')
}  # past
par(mfrow=c(1,1))

# Present Model #
minPred.present <- minValue(present.study.region)
maxPred.present <- maxValue(present.study.region)
names(minPred.present) <- names(maxPred.present) <- names(present.study.region)

medianPred.present <- apply(records.present[ , predictors], 2, median)
medianPred.present

env.present <- as.data.frame(medianPred.present)
env.present <- t(env.present)
env.present <- env.present[rep(1, 100), ]
row.names(env.present) <- 1:nrow(env.present)
head(env.present)

par(mfrow=c(2,3), mar=c(4,4,1,1), ps=8, family="sans")
for (pred in predictors) {
  # make copy of data frame
  thisEnv <- env.present
  # now vary focal predictor from min to max value... 
  # all other predictors keep median value
  thisEnv[ , pred] <- seq(minPred.present[pred], maxPred.present[pred], length.out=100)
  # make prediction using this data frame
  prediction <- predict(Model.raw.present, thisEnv, type='cloglog')
  # plot
  plot(x=thisEnv[ , pred], y=prediction, ylim=c(0, 1), xlab=pred,
       ylab='Suitability', main="", type='l', col='black',
       lty=1, lwd=1, cex.lab=1.3)
  # add species' presences (top rug)
  rug(records.present[ , pred], side=3, col='darkblue')
  # add background sites (bottom rug)
  rug(targetBg.raw.present[ , pred], side=1, col='darkblue')
}  # present, same features and regularization as past model (lpq 0.5)


## Predictions and maps ##

# Project models
tunedMap.raw.past <- predict(past.study.region[[predictors]], tunedModel.raw.past, type='cloglog') # Past Model on past climate
tunedMap.raw.present <- predict(present.study.region[[predictors]], tunedModel.raw.past, type='cloglog') # Past Model on present climate
Map.raw.now <- predict(present.study.region[[predictors]], Model.raw.present, type='cloglog') # Present Model on present climate

### FIGURE 2 (note that final map has been reprojected in QGIS)

#pdf("Results/Figure2.pdf", width=6.614173, height=4.55, pointsize=8) # to save the figure as pdf; page width 168mm
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))

plot(crop(tunedMap.raw.past,extEU), main='', 
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous")) # Past Model projected on past climate, occurrences in historic native range
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.past, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.1)
plot(crop(tunedMap.raw.present,extEU), main="", 
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F) # Past Model projected on present climate, occurrences in historic native range
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.past, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.1)
plot(crop(tunedMap.raw.present,extEU), main="", 
       breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F) # same, but with all occurrences
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.present, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.1)
plot(crop(Map.raw.now,extEU), main="", 
       breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F) # Present Model projected on present climate, all occurrences
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.present, pch=21, bg=alpha('mediumseagreen', 0.5), cex=.1)
  
dev.off()

## FIGURE S2: entire study area 
#pdf("Results/FigureS2.pdf", width=8, height=5.5)
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))

plot(tunedMap.raw.past, main='',
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous")) # Past Model projected on past climate, occurrences in historic native range
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.past, pch=20, cex=.005)
plot(tunedMap.raw.present, main="", 
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F) # Past Model projected on present climate, occurrences in historic native range
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.past, pch=20, cex=.005)
plot(tunedMap.raw.present, main="", 
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F)  # same, but with all occurrences
  lines(native.range, col="grey20")
  points(latitude ~ longitude, records.present, pch=20, cex=.005)
plot(Map.raw.now, main="", 
     breaks=seq(0,1,.01), col=cal_palette("figmtn", n=100, type="continuous"), legend=F) # Present Model projected on present climate, all occurrences
lines(native.range, col="grey20")
points(latitude ~ longitude, records.present, pch=20, cex=.005)

dev.off()

# Export map elements for map formatting in QGIS
writeRaster(tunedMap.raw.past, "Results/finalmap_past_past.tif", format="GTiff")
writeRaster(tunedMap.raw.present, "Results/finalmap_past_present.tif", format="GTiff")
writeRaster(Map.raw.now, "Results/finalmap_present_present.tif", format="GTiff")

write.csv(records.past, "Results/records_past.csv", row.names=F)
write.csv(records.present, "Results/records_present.csv", row.names=F)
write.table(as.character(cal_palette("figmtn", n=100, type="continuous")), "Results/calpalette_figmtn.txt", col.names=F, row.names=F) # colors


#### ====================== Thresholded maps (Appendix S7) ======================================

# Past model to past climate
predPres.pastmodel.past <- extract(tunedMap.raw.past, cbind(records.past$longitude, records.past$latitude))
predBg.pastmodel.past <- extract(tunedMap.raw.past, cbind(targetBg.raw.past$longitude, targetBg.raw.past$latitude))

# Past model to present climate
predPres.pastmodel <- extract(tunedMap.raw.present, cbind(records.present$longitude, records.present$latitude))
predBg.pastmodel <- extract(tunedMap.raw.present, cbind(targetBg.raw.present$longitude, targetBg.raw.present$latitude))

# Present model to past climate (for sensitivity only)
Map.raw.now.past <- predict(present.study.region[[predictors]], Model.raw.present, type='cloglog') # Present Model on past climate
predPres.presentmodel.past <- extract(Map.raw.now.past, cbind(records.past$longitude, records.past$latitude))
predBg.presentmodel.past <- extract(Map.raw.now.past, cbind(targetBg.raw.past$longitude, targetBg.raw.past$latitude))

# Present model to present climate
predPres.presentmodel <- extract(Map.raw.now, cbind(records.present$longitude, records.present$latitude))
predBg.presentmodel <- extract(Map.raw.now, cbind(targetBg.raw.present$longitude, targetBg.raw.present$latitude))

# create evaluation object for each model
eval.pastmodel.past <- evaluate(p=as.vector(predPres.pastmodel.past), a=as.vector(predBg.pastmodel.past), tr=seq(0, 1, by=0.01))
eval.pastmodel.past
eval.pastmodel <- evaluate(p=as.vector(predPres.pastmodel), a=as.vector(predBg.pastmodel), tr=seq(0, 1, by=0.01))
eval.pastmodel
eval.presentmodel <- evaluate(p=as.vector(predPres.presentmodel), a=as.vector(predBg.presentmodel), tr=seq(0, 1, by=0.01))
eval.presentmodel

# threshold with maximum across past and present model
threshold.maxboth <- eval.pastmodel@t[which.max(eval.pastmodel@TPR + eval.pastmodel@TNR + eval.presentmodel@TPR + eval.presentmodel@TNR)]
threshold.maxboth # 0.37 
threshold <- threshold.maxboth

threshold.past.past <- tunedMap.raw.past >= threshold # past model in the past  
threshold.past <- tunedMap.raw.present >= threshold # past model in the present
threshold.present <- Map.raw.now >= threshold # present threshold

# sensitivity 
sum(predPres.pastmodel.past >= threshold) / length(predPres.pastmodel.past) # sensitivity of past model in the past
sum(predPres.pastmodel >= threshold) / length(predPres.pastmodel) # sensitivity of past model in the present
sum(predPres.presentmodel.past >= threshold) / length(predPres.presentmodel.past) # sensitivity of present model in the past
sum(predPres.presentmodel >= threshold) / length(predPres.presentmodel) # sensitivity of present model in the present

# plot thresholded maps by themselves
plot(crop(threshold.past, extEU))
  lines(native.range)
  points(latitude ~ longitude, records.past, pch=16, cex=.5) # historic native range (present climate)
plot(crop(threshold.present, extEU))
  lines(native.range)
  points(latitude ~ longitude, records.present, pch=16, cex=.5) # expanded native range (present climate)

# combine thresholded maps
threshold.past[threshold.past>0] <- -1
threshold.present[threshold.present>0] <- 2
threshold.difference <- threshold.past + threshold.present
# unsuitable habitat has the value 0
# outcome of adding up the past and present rasters, with suitable habitat in:
  # neither model:      0 + 0 = 0 
  # past model only:   -1 + 0 = -1  
  # both models:       -1 + 2 = 1  
  # present model only: 0 + 2 = 2  

### FIGURE S5 ###
#pdf("Results/FigureS5.pdf", width=6.614173, height=4.55, pointsize=8)
plot(crop(threshold.difference, extEU), 
     col=cal_palette("figmtn", n=5, type="discrete"))
 lines(study.region)
 lines(native.range, col="grey30", lwd=2) 
 points(latitude ~ longitude, records.present, pch=16, cex=.5, col="white") # expanded native range
 points(latitude ~ longitude, records.past, pch=16, cex=.5) # historic native range
dev.off()

# Save to QGIS for map formatting
writeRaster(threshold.difference, "Results/thresholdFigS5.tif", format="GTiff")


#### ====================== Model evaluation ======================================

### Random cross-validation (k-folds): PAST ###

# calculate k-folds for presences and background sites
kPres.past <- kfold(x=records.past, k=5)
kBg.past <- kfold(x=targetBg.raw.past, k=5)
head(kPres.past)
head(kBg.past)

length(kPres.past) - table(kPres.past) # training presences in each k-fold (in the model)
table(kPres.past) # test presences in each k-fold

length(kBg.past) - table(kBg.past) # training background points in each k-fold
table(kBg.past) # test background points

# map
plot(crop(past.study.region, extEU)[[1]], main='k-fold #1')
points(records.past$longitude, records.past$latitude)
points(records.past$longitude[kPres.past==1],
       records.past$latitude[kPres.past==1],
       bg='red',
       pch=21)
legend('topleft',
       legend=c('Training presence', 'Test presence'),
       pch=c(1, 16),
       col=c('black', 'red'),
       bg='white',
       cex=0.8)

# for storing AUC and CBI
aucRandom.past <- cbiRandom.past <- rep(NA, 5)

# cycle through each k-fold
for (i in 1:5) {
  
  say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    records.past[kPres.past!=i, predictors],
    targetBg.raw.past[kBg.past!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(kPres.past!=i)), rep(0, sum(kBg.past!=i)))
  
  trainData <- cbind(presBg, envData)
  
  # tuned model
  model <- trainMaxNet(
    data=trainData,
    regMult=0.5,
    classes='lpq',
    testClasses=F,
    verbose=FALSE
  )
  
  # predict to presences and background sites
  predPres <- raster::predict(model, newdata=records.past[kPres.past==i, ], type='cloglog')
  predBg <- raster::predict(model, newdata=targetBg.raw.past[kBg.past==i, ], type='cloglog')
  
  # evaluate and remember result
  thisEval <- dismo::evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisAuc <- thisEval@auc
  thisCbi <- contBoyce(pres=predPres, bg=predBg)
  
  say(': AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2))
  
  aucRandom.past[i] <- thisAuc
  cbiRandom.past[i] <- thisCbi
  
}

# Results for each k-fold, overall mean and standard deviation
aucRandom.past
mean(aucRandom.past); sd(aucRandom.past)

cbiRandom.past
mean(cbiRandom.past); sd(cbiRandom.past)


### Random cross-validation (k-folds): PAST MODEL ON PRESENT DATA ###

# for storing AUC and CBI
aucRandom.pasttopresent <- cbiRandom.pasttopresent <- rep(NA, 5)

# cycle through each k-fold
for (i in 1:5) {
  
  say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    records.past[kPres.past!=i, predictors],
    targetBg.raw.past[kBg.past!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(kPres.past!=i)), rep(0, sum(kBg.past!=i)))
  
  trainData <- cbind(presBg, envData)
  
  # tuned model
  model <- trainMaxNet(
    data=trainData,
    regMult=0.5,
    classes='lpq',
    testClasses=F,
    verbose=FALSE
  )
  
  # predict to presences and background sites of the PRESENT 
  predPres <- raster::predict(model, newdata=records.present[kPres.present==i, ], type='cloglog')
  predBg <- raster::predict(model, newdata=targetBg.raw.present[kBg.present==i, ], type='cloglog')
  
  # evaluate and remember result
  thisEval <- dismo::evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisAuc <- thisEval@auc
  thisCbi <- contBoyce(pres=predPres, bg=predBg)
  
  say(': AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2))
  
  aucRandom.pasttopresent[i] <- thisAuc
  cbiRandom.pasttopresent[i] <- thisCbi
  
}

# Results for each k-fold, overall mean and standard deviation
aucRandom.pasttopresent
mean(aucRandom.pasttopresent); sd(aucRandom.pasttopresent)

cbiRandom.pasttopresent
mean(cbiRandom.pasttopresent); sd(cbiRandom.pasttopresent)


### Random cross-validation (k-folds): PRESENT MODEL ###

# calculate k-folds for presences and background sites
kPres.present <- kfold(x=records.present, k=5)
kBg.present <- kfold(x=targetBg.raw.present, k=5)
head(kPres.present)
head(kBg.present)

length(kPres.present) - table(kPres.present) # training presences in each k-fold (in the model)
table(kPres.present) # test presences in each k-fold

length(kBg.present) - table(kBg.present) # training background points in each k-fold
table(kBg.present) # test background points

# map
plot(crop(past.study.region, extEU)[[1]], main='k-fold #1')
points(records.present$longitude, records.present$latitude)
points(records.present$longitude[kPres.present==1],
       records.present$latitude[kPres.present==1],
       bg='red',
       pch=21)
legend('topleft',
       legend=c('Training presence', 'Test presence'),
       pch=c(1, 16),
       col=c('black', 'red'),
       bg='white',
       cex=0.8)

# for storing AUC and CBI
aucRandom.present <- cbiRandom.present <- rep(NA, 5)

# cycle through each k-fold
for (i in 1:5) {
  
  say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    records.present[kPres.present!=i, predictors],
    targetBg.raw.present[kBg.present!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(kPres.present!=i)), rep(0, sum(kBg.present!=i)))
  
  trainData <- cbind(presBg, envData)
  
  # tuned model
  model <- trainMaxNet(
    data=trainData,
    regMult=0.5, # only use 1 regularization parameter (the same one as the main models)
    classes='lpq',
    testClasses=F, # don't try out different features, keep lpq fixed
    verbose=FALSE
  )
  
  # predict to presences and background sites
  predPres <- raster::predict(model, newdata=records.present[kPres.present==i, ], type='cloglog')
  predBg <- raster::predict(model, newdata=targetBg.raw.present[kBg.present==i, ], type='cloglog')
  
  # evaluate and remember result
  thisEval <- dismo::evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisAuc <- thisEval@auc
  thisCbi <- contBoyce(pres=predPres, bg=predBg)
  
  say(': AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2))
  
  aucRandom.present[i] <- thisAuc
  cbiRandom.present[i] <- thisCbi
  
}

# Results for each k-fold, overall mean and standard deviation
aucRandom.present
mean(aucRandom.present); sd(aucRandom.present)

cbiRandom.present
mean(cbiRandom.present); sd(cbiRandom.present)

### Store results in Table 3 ###
modelevaluation <- data.frame("model" = c("past", "past", "present"),
                              "project.to" = c("past", "present", "present"),
                              "n.training.pres" = rep(NA,3),
                              "n.test.pres" = rep(NA,3),
                              "n.training.bg" = rep(NA,3),
                              "n.test.bg" = rep(NA,3),
                              "AUC mean" = rep(NA,3),
                              "AUC sd" = rep(NA,3),
                              "CBI mean" = rep(NA,3),
                              "CBI sd" = rep(NA,3))
modelevaluation

modelevaluation[1,3:4] <- c((length(kPres.past) - table(kPres.past))[1], table(kPres.past)[1])
modelevaluation[1,5:6] <- c((length(kBg.past) - table(kBg.past))[1], table(kBg.past)[1])
modelevaluation[2,3:4] <- c((length(kPres.past) - table(kPres.past))[1], table(kPres.present)[1])
modelevaluation[2,5:6] <- c((length(kBg.past) - table(kBg.past))[1], table(kBg.present)[1])
modelevaluation[3,3:4] <- c((length(kPres.present) - table(kPres.present))[1], table(kPres.present)[1])
modelevaluation[3,5:6] <- c((length(kBg.present) - table(kBg.present))[1], table(kBg.present)[1])

modelevaluation[1,7:8] <- c(mean(aucRandom.past), sd(aucRandom.past))
modelevaluation[2,7:8] <- c(mean(aucRandom.pasttopresent), sd(aucRandom.pasttopresent))
modelevaluation[3,7:8] <- c(mean(aucRandom.present), sd(aucRandom.present))

modelevaluation[1,9:10] <- c(mean(cbiRandom.past), sd(cbiRandom.past))
modelevaluation[2,9:10] <- c(mean(cbiRandom.pasttopresent), sd(cbiRandom.pasttopresent))
modelevaluation[3,9:10] <- c(mean(cbiRandom.present), sd(cbiRandom.present))

modelevaluation[,7:10] <- round(modelevaluation[,7:10],2)

modelevaluation # Table 3; note that results may vary slightly due to random k-folds.


#### ====================== Difference between projections in the invaded range =====================

### California ###

# environmental predictors plotted in California
plot(present.bioclim.CA[[predictors]])

# model projections
CAmap.raw.past <- predict(present.bioclim.CA[[predictors]], tunedModel.raw.past, type='cloglog') # Past Model on present climate
CAmap.raw.now <- predict(present.bioclim.CA[[predictors]], Model.raw.present, type='cloglog') # Present Model on present climate

# Minimum presence threshold: 
# find all other areas that are at least as suitable, as the least-suitable habitat already invaded in California
predPres.CA.past <- extract(CAmap.raw.past, cbind(calflora.records$x, calflora.records$y))
predPres.CA.now <- extract(CAmap.raw.now, cbind(calflora.records$x, calflora.records$y))

CAmap.raw.past.minPres <- CAmap.raw.past >= min(predPres.CA.past, na.rm=T)  # past threshold
CAmap.raw.now.minPres <- CAmap.raw.now >= min(predPres.CA.now, na.rm=T) # present threshold

# create map
CAmap.raw.past.minPres[CAmap.raw.past.minPres>0] <- -1 # suitable habitat in past model gets the value -1
CAmap.raw.now.minPres[CAmap.raw.now.minPres>0] <- 2    # suitable habitat in present model gets the value 2
CAmap.difference.minPres <- CAmap.raw.past.minPres + CAmap.raw.now.minPres # add up both maps
  # unsuitable habitat has the value 0
  # outcome of adding up the past and present rasters, with suitable habitat in:
    # neither model:      0 + 0 = 0 
    # past model only:   -1 + 0 = -1  
    # both models:       -1 + 2 = 1  
    # present model only: 0 + 2 = 2  

# compute Multivariate Environmental Suitability Surfaces (MESS) to find nonanalog climates
mess.CA.past <- mess(present.bioclim.CA[[predictors]], trainData.raw.past)
mess.CA.present <- mess(present.bioclim.CA[[predictors]], trainData.raw.present)

plot(mess.CA.present<0) # cells with MESS<0 will be considered nonanalog

CAmap.difference.minPres[mess.CA.past<0] <- -2 # assign the value -2 to cells outside mess
CAmap.difference.minPres[mess.CA.present<0] <- -2

# Figure 3c
plot(CAmap.difference.minPres, 
     col=c("grey95", "#7bb1c7", "#cceefc", "#8ac7e0", "#9adef9"), # outside mess, past, neither, both, present
     legend=F, main="")
points(y ~ x, calflora.records, pch=16)
lines(USstates)

### Australia ###

# environmental predictors plotted in Australia
plot(present.bioclim.AU[[predictors]])

# model projections
AUmap.raw.past <- predict(present.bioclim.AU[[predictors]], tunedModel.raw.past, type='cloglog') # past model on present climate
AUmap.raw.now <- predict(present.bioclim.AU[[predictors]], Model.raw.present, type='cloglog')    # present model on present climate

# Minimum presence treshold 
predPres.AU.past <- extract(AUmap.raw.past, cbind(Australia.records$x, Australia.records$y))
predPres.AU.now <- extract(AUmap.raw.now, cbind(Australia.records$x, Australia.records$y))

AUmap.raw.past.minPres <- AUmap.raw.past >= min(predPres.AU.past, na.rm=T) # past threshold
AUmap.raw.now.minPres <- AUmap.raw.now >= min(predPres.AU.now, na.rm=T) # present threshold

# create map
AUmap.raw.past.minPres[AUmap.raw.past.minPres>0] <- -1
AUmap.raw.now.minPres[AUmap.raw.now.minPres>0] <- 2
AUmap.difference.minPres <- AUmap.raw.past.minPres + AUmap.raw.now.minPres

# compute MESS surfaces
mess.AU.past <- mess(present.bioclim.AU[[predictors]], trainData.raw.past) # past training data, correspondence to present climate of Australia
mess.AU.present <- mess(present.bioclim.AU[[predictors]], trainData.raw.present) # present training data, correspondence to present climate of Australia

plot(mess.AU.past<0)  
  lines(Australia)

AUmap.difference.minPres[mess.AU.past<0] <- -2 # assign the value -2 to cells outside mess
AUmap.difference.minPres[mess.AU.present<0] <- -2

# Figure 3f
plot(crop(AUmap.difference.minPres, extAU), 
     col=c("grey95", "#7bb1c7", "#cceefc", "#8ac7e0", "#9adef9"), # outside mess, past, neither, both, present
     legend=F, main="")
points(y ~ x, Australia.records,cex=.4, pch=16) # black outline for points
lines(Australia, lwd=1.5)

