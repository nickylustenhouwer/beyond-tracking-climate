# R code for:
# Beyond tracking climate: niche evolution during native range expansion and its implications for novel invasions
# N. Lustenhouwer & I.M. Parker

# Script 3: COUE analysis

# Last edit: August 27, 2021

# Code adapted from:
  # Di Cola, V. et al. 2017. ecospat: an R package to support spatial analyses and modeling of species niches and distributions. - Ecography 40: 774–787


## Load climate and occurrence data from starting script
source("Rscripts/LustenhouwerParker_start.R")

## Parameter settings used for all analyses
bioclimsub <- predictors.numeric + 2 # column numbers of the environmental predictors (columns 1 and 2 are x and y coordinates)
cutoff <- 0.1 # intersection cutoff: the quantile of the environmental density used to remove marginal climates


### ================== MODIFICATIONS TO THE PLOTTING FUNCTIONS FROM THE ECOSPAT PACKAGE ==========

# a few minor modifications to plot formatting (color etc)
# no changes to the actual plot

## edits to ecospat.plot.niche.dyn {ecospat} ##
# two edits (indicated below) to change the colors and transparancy

myplot <- function (z1, z2, quant, title = "", name.axis1 = "Axis 1", name.axis2 = "Axis 2", 
                    interest = 1, colz1 = "#00FF0050", colz2 = "#FF000050", colinter = "#0000FF50", 
                    colZ1 = "green3", colZ2 = "red3") 
{
  if (is.null(z1$y)) {
    R <- length(z1$x)
    x <- z1$x
    xx <- sort(rep(1:length(x), 2))
    y1 <- z1$z.uncor/max(z1$z.uncor)
    Y1 <- z1$Z/max(z1$Z)
    if (quant > 0) {
      Y1.quant <- quantile(z1$Z[which(z1$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z1$Z)
    }
    else {
      Y1.quant <- 0
    }
    Y1.quant <- Y1 - Y1.quant
    Y1.quant[Y1.quant < 0] <- 0
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 
                                           2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 
                                           2)]
    y2 <- z2$z.uncor/max(z2$z.uncor)
    Y2 <- z2$Z/max(z2$Z)
    if (quant > 0) {
      Y2.quant <- quantile(z2$Z[which(z2$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z2$Z)
    }
    else {
      Y2.quant = 0
    }
    Y2.quant <- Y2 - Y2.quant
    Y2.quant[Y2.quant < 0] <- 0
    yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2) * 
                                           2)]
    YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2) * 
                                           2)]
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = colz1, border = 0)
    polygon(x[xx], c(0, y2[yy2], 0, 0), col = colz2, border = 0)
    polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, 
                              min, na.exclude = TRUE), 0, 0), col = colinter, border = 0)
    lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, 
          lty = "dashed")
    lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, 
          lty = "dashed")
    lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
    lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
    segments(x0 = 0, y0 = 0, x1 = max(x[xx]), y1 = 0, col = "white")
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, col = "white")
    seg.cat <- function(inter, cat, col.unf, col.exp, col.stab) {
      if (inter[3] == 0) {
        my.col = 0
      }
      if (inter[3] == 1) {
        my.col = col.unf
      }
      if (inter[3] == 2) {
        my.col = col.stab
      }
      if (inter[3] == -1) {
        my.col = col.exp
      }
      segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2], 
               col = my.col, lwd = 4, lty = 2)
    }
    cat <- ecospat.niche.dyn.index(z1, z2, intersection = quant)$dyn
    inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1])
    apply(inter, 1, seg.cat, col.unf = "#00FF0050", col.exp = "#FF000050", 
          col.stab = "#0000FF50")
  }
  if (!is.null(z1$y)) {
    z <- t(as.matrix(z1$w + 2 * z2$w))[, nrow(as.matrix(z1$z.uncor)):1]
    z1$Z <- t(as.matrix(z1$Z))[, nrow(as.matrix(z1$Z)):1]
    z2$Z <- t(as.matrix(z2$Z))[, nrow(as.matrix(z2$Z)):1]
    if (interest == 1) {
      image(x = z1$x, y = z1$y, z = z, col = c("#FFFFFF00", colz1, colz2, colinter), # edit 1: switched order of this image and the one below
            xlab = name.axis1, ylab = name.axis2)
      image(x = z1$x, y = z1$y, z = t(as.matrix(z1$z.uncor))[, 
            nrow(as.matrix(z1$z.uncor)):1], col = paste(gray(100:0/100),"99", sep=""), # edit 2: added "99" to hex color codes for transparancy
            zlim = c(1e-05, cellStats(z1$z.uncor, "max")), 
            add=T)
    }
    if (interest == 2) {
      image(x = z2$x, y = z2$y, z = z, col = c("#FFFFFF00", colz1, colz2, colinter), # same as edit 1
            xlab = name.axis1, ylab = name.axis2)
      image(x = z2$x, y = z2$y, z = t(as.matrix(z2$z.uncor))[, 
            nrow(as.matrix(z2$z.uncor)):1], col = paste(gray(100:0/100),"99",sep=""), # same as edit 2
            zlim = c(1e-05, cellStats(z2$z.uncor, "max")), 
            add = T)
    }
    title(title)
    contour(x = z1$x, y = z1$y, z1$Z, add = TRUE, levels = quantile(z1$Z[z1$Z > 
            0], c(0, quant)), drawlabels = FALSE, lty = c(1,2), col = colZ1)
    contour(x = z2$x, y = z2$y, z2$Z, add = TRUE, levels = quantile(z2$Z[z2$Z > 0], 
            c(0, quant)), drawlabels = FALSE, lty = c(1,2), col = colZ2)
  }
}


## edits to ecospat.shift.centroids {ecospat} ##
# changed the linetype and selected individuals colors for the climate and species arrow

my.centroid <- function (sp1, sp2, clim1, clim2, col.sp = "red", col.clim = "black") # added 'col.sp' and 'col.clim' instead of one 'col' argument
{
  if (ncol(as.matrix(sp1)) == 2) {
    arrows(median(sp1[, 1]), median(sp1[, 2]), median(sp2[, 
        1]), median(sp2[, 2]), col = col.sp, lwd = 2, length = 0.1) # col = col.sp
    arrows(median(clim1[, 1]), median(clim1[, 2]), median(clim2[, 
        1]), median(clim2[, 2]), lty = 1, col = col.clim, lwd = 2, # lty= 1 instead of 11, col = col.clim
           length = 0.1)
  }
  else {
    arrows(median(sp1), 0.025, median(sp2), 0.025, col = col.sp, # col = col.sp 
           lwd = 2, length = 0.1)
    arrows(median(clim1), -0.025, median(clim2), -0.025, 
           lty=1, col = col.clim, lwd = 2, length = 0.1) # lty = 1 instead of 11, col = col.clim
  }
}



#### ================= COMPARING TIME PERIODS IN THE NATIVE RANGE ================== 

### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both time periods
pca.env.pastVpresent <- dudi.pca(rbind(past.occ.original, present.occ)[,bioclimsub], scannf=FALSE, nf=2)  #  this calibrates the PCA in the entire study area (past and present)
ecospat.plot.contrib(contrib=pca.env.pastVpresent$co, eigen=pca.env.pastVpresent$eig) # plot variables contribution

# predict the scores on the axes
scores.globclim.pastVpresent <- pca.env.pastVpresent$li # PCA scores for the whole study area 
scores.sp.past.pastVpresent <- suprow(pca.env.pastVpresent, past.occ.original[which(past.occ.original[,"species_occ"]==1),bioclimsub])$li # PCA scores for the past 
head(scores.sp.past.pastVpresent)
scores.sp.present.pastVpresent <- suprow(pca.env.pastVpresent, present.occ[which(present.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for the present 
scores.clim.past.pastVpresent <- suprow(pca.env.pastVpresent, past.occ.original[,bioclimsub])$li # PCA scores for the whole past study area 
scores.clim.present.pastVpresent <- suprow(pca.env.pastVpresent, present.occ[,bioclimsub])$li # PCA scores for the whole present study area

## Calculate the Occurrence Densities Grid (in environmental space, 100 x 100 cells)

# For the past
grid.clim.past.pastVpresent <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVpresent,
                                                     glob1=scores.clim.past.pastVpresent, sp=scores.sp.past.pastVpresent, R=100, th.sp=0)

# For the present
grid.clim.present.pastVpresent <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVpresent,
                                                        glob1=scores.clim.present.pastVpresent, sp=scores.sp.present.pastVpresent, R=100, th.sp=0)  

## Calculate Niche Overlap (D value) 
D.pastVpresent <- ecospat.niche.overlap(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent, cor=TRUE)$D
D.pastVpresent

## Niche Equivalency Test
eq.test.pastVpresent <- ecospat.niche.equivalency.test(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent,
                                                       rep=100, alternative = "lower") 

ecospat.plot.overlap.test(eq.test.pastVpresent, "D", "Equivalency")
# H0 = niche equivalency (niche overlap is constant when shifting occurrences around between ranges)
# alternative=lower tests for niche divergence (niches less similar than random: D is lower than random)
# alternative=greater tests for niche conservatism (niches more similar than random: D is larger than random)

# Result: niches became different during range expansion (D is lower than the null distribution).
eq.test.pastVpresent$p.D

## Niche Similarity Test
# shifting both niches (rand.type=1)
sim.test.pastVpresent <- ecospat.niche.similarity.test(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent,
                                                       rep=100, alternative = "greater", rand.type = 1)  

ecospat.plot.overlap.test(sim.test.pastVpresent, "D", "Similarity")
# H0 = niche similarity (niches are as similar to each other as expected by chance)
# To test for niche conservatism, use the alternative ‘greater’, i.e. the niche overlap is more similar than random expectations.

# Result: niches are still more similar to each other than expected by chance (compared to any random niche from the range)
sim.test.pastVpresent$p.D

## Delimiting niche categories and quantifying niche dynamics in analogue climates (expansion/stability/unfilling)
w.pastVpresent <- ecospat.niche.dyn.index(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent, intersection=cutoff)$dynamic.index.w
w.pastVpresent

# Visualizing niche categories, niche dynamics and climate analogy between time periods
myplot(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent, quant=cutoff, 
       interest=1, # interest=1 plots past, interest=2 plots present density
       title= "Niche Overlap", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988",
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
my.centroid(scores.sp.past.pastVpresent, scores.sp.present.pastVpresent, scores.clim.past.pastVpresent, scores.clim.present.pastVpresent, 
            col.sp="white", col.clim="black") # shift of niche centroid and climate centroid
# purple shading = niche space occupied in both time periods (stability)
# orange shading = past only (unfilling)
# blue shading = present only (expansion)
# dark shading = density of occurrences
# orange and blue lines show percentiles of the global climate (Eurasian Holarctic) in the past and present, respectively (solid=100%, dashed=90%)

# Plot Variables Contribution
ecospat.plot.contrib(contrib=pca.env.pastVpresent$co, eigen=pca.env.pastVpresent$eig)

## Figure 1 ##

#pdf("Figure1.pdf", width=8.3, height=4.7)
layout(matrix(c(1,1,1,1,2,3,2,4), nrow=2, ncol=4))
myplot(grid.clim.past.pastVpresent, grid.clim.present.pastVpresent, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/blue
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
my.centroid(scores.sp.past.pastVpresent, scores.sp.present.pastVpresent, scores.clim.past.pastVpresent, scores.clim.present.pastVpresent, 
            col.sp="white", col.clim="black") # shift of niche centroid and climate centroid
ecospat.plot.contrib(contrib=pca.env.pastVpresent$co, eigen=pca.env.pastVpresent$eig)
ecospat.plot.overlap.test(eq.test.pastVpresent, "D", "Equivalency")
ecospat.plot.overlap.test(sim.test.pastVpresent, "D", "Similarity")
dev.off()

# Edits in inkscape: increased label sizes on PCA plot, added panel labels and axis labels, removed unnecessary text (panel headers and p.values).


#### ================ PAST NATIVE RANGE CLIMATE AND OCCURRENCES TO CALIFORNIA ===================

### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both Eurasia and North America
pca.env.pastVca <- dudi.pca(rbind(past.occ.original, calflora.occ)[,bioclimsub], scannf=FALSE, nf=2)  #  this calibrates the PCA in the entire study area (Eurasian and North American Holarctic)
ecospat.plot.contrib(contrib=pca.env.pastVca$co, eigen=pca.env.pastVca$eig) # plot variables contribution

# predict the scores on the axes
scores.globclim.pastVca <- pca.env.pastVca$li # PCA scores for the whole study area 
scores.sp.past.pastVca <- suprow(pca.env.pastVca, past.occ.original[which(past.occ.original[,"species_occ"]==1),bioclimsub])$li # PCA scores for the past, native distribution 
head(scores.sp.past.pastVca)
scores.sp.present.pastVca <- suprow(pca.env.pastVca, calflora.occ[which(calflora.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for California
scores.clim.past.pastVca <- suprow(pca.env.pastVca, past.occ.original[,bioclimsub])$li # PCA scores for all of Eurasia (past climate)
scores.clim.present.pastVca <- suprow(pca.env.pastVca, calflora.occ[,bioclimsub])$li # PCA scores for all of North America

# Calculate the Occurrence Densities Grid (in environmental space, 100 x 100 cells)

# For the native range (past model)
grid.clim.past.pastVca <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVca,
                                                glob1=scores.clim.past.pastVca, sp=scores.sp.past.pastVca, R=100, th.sp=0)

# For North America
grid.clim.present.pastVca <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVca,
                                                   glob1=scores.clim.present.pastVca, sp=scores.sp.present.pastVca, R=100, th.sp=0)  

# Calculate Niche Overlap (D value)
D.pastVca <- ecospat.niche.overlap (grid.clim.past.pastVca, grid.clim.present.pastVca, cor=TRUE)$D
D.pastVca

## Niche Equivalency Test
eq.test.pastVca <- ecospat.niche.equivalency.test(grid.clim.past.pastVca, grid.clim.present.pastVca,
                                                  rep=100, alternative = "lower") 

ecospat.plot.overlap.test(eq.test.pastVca, "D", "Equivalency")
eq.test.pastVca$p.D 

## Niche Similarity Test
# rand.type=2: shifting randomly the invasive niche in the invaded study area
sim.test.pastVca <- ecospat.niche.similarity.test(grid.clim.past.pastVca, grid.clim.present.pastVca,
                                                  rep=100, alternative = "greater", rand.type = 2)  

ecospat.plot.overlap.test(sim.test.pastVca, "D", "Similarity")
sim.test.pastVca$p.D 
# Result: niches are not equivalent, but we cannot rule out that the niche shift is due to differences in available habitat.

## Delimiting niche categories and quantifying niche dynamics in analogue climates (expansion/stability/unfilling)
w.pastVca <- ecospat.niche.dyn.index(grid.clim.past.pastVca, grid.clim.present.pastVca, intersection=cutoff)$dynamic.index.w
w.pastVca

# Visualizing niche categories, niche dynamics and climate analogy between ranges
myplot(grid.clim.past.pastVca, grid.clim.present.pastVca, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "Niche Overlap", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/green
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
# Plot Variables Contribution
ecospat.plot.contrib(contrib=pca.env.pastVca$co, eigen=pca.env.pastVca$eig)


#### ================ PRESENT NATIVE RANGE CLIMATE AND OCCURRENCES TO CALIFORNIA ===================

### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both Eurasia and North America
pca.env.presentVca <- dudi.pca(rbind(present.occ, calflora.occ)[,bioclimsub], scannf=FALSE, nf=2)  #  this calibrates the PCA in the entire study area (Eurasian and North American Holarctic)
ecospat.plot.contrib(contrib=pca.env.presentVca$co, eigen=pca.env.presentVca$eig) # plot variables contribution

# predict the scores on the axes
scores.globclim.presentVca <- pca.env.presentVca$li # PCA scores for the whole study area 
scores.sp.past.presentVca <- suprow(pca.env.presentVca, present.occ[which(present.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for the present native distribution
head(scores.sp.past.presentVca)
scores.sp.present.presentVca <- suprow(pca.env.presentVca, calflora.occ[which(calflora.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for occurrences in California
scores.clim.past.presentVca <- suprow(pca.env.presentVca, present.occ[,bioclimsub])$li # PCA scores for the whole native range, present climate (Eurasian Holarctic)
scores.clim.present.presentVca <- suprow(pca.env.presentVca, calflora.occ[,bioclimsub])$li # PCA scores for all of North America

# Calculate the Occurrence Densities Grid (in environmental space, 100 x 100 cells)

# For the native range (present climate)
grid.clim.past.presentVca <- ecospat.grid.clim.dyn(glob=scores.globclim.presentVca,
                                                   glob1=scores.clim.past.presentVca, sp=scores.sp.past.presentVca, R=100, th.sp=0)

# For North America
grid.clim.present.presentVca <- ecospat.grid.clim.dyn(glob=scores.globclim.presentVca,
                                                      glob1=scores.clim.present.presentVca, sp=scores.sp.present.presentVca, R=100, th.sp=0)  

# Calculate Niche Overlap (D value)
D.presentVca <- ecospat.niche.overlap (grid.clim.past.presentVca, grid.clim.present.presentVca, cor=TRUE)$D
D.presentVca

## Niche Equivalency Test
eq.test.presentVca <- ecospat.niche.equivalency.test(grid.clim.past.presentVca, grid.clim.present.presentVca,
                                                     rep=100, alternative = "lower") 

ecospat.plot.overlap.test(eq.test.presentVca, "D", "Equivalency")
eq.test.presentVca$p.D

## Niche Similarity Test
sim.test.presentVca <- ecospat.niche.similarity.test(grid.clim.past.presentVca, grid.clim.present.presentVca,
                                                     rep=100, alternative = "greater", rand.type = 2)  

ecospat.plot.overlap.test(sim.test.presentVca, "D", "Similarity")
sim.test.presentVca$p.D

# Delimiting niche categories and quantifying niche dynamics in analogue climates (expansion/stability/unfilling)
w.presentVca <- ecospat.niche.dyn.index(grid.clim.past.presentVca, grid.clim.present.presentVca, intersection=cutoff)$dynamic.index.w
w.presentVca

# Visualizing niche categories, niche dynamics and climate analogy between ranges
myplot(grid.clim.past.presentVca, grid.clim.present.presentVca, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "Niche Overlap", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/green
       colZ1 = "#e29f60", colZ2 = "#52c6f5")

# Plot Variables Contribution
ecospat.plot.contrib(contrib=pca.env.presentVca$co, eigen=pca.env.presentVca$eig)


#### ================ PAST NATIVE RANGE CLIMATE AND OCCURRENCES TO AUSTRALIA ===================

### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both Eurasia and Australia
pca.env.pastVau <- dudi.pca(rbind(past.occ.original, Australia.occ)[,bioclimsub], scannf=FALSE, nf=2)  #  this calibrates the PCA in the entire study area (Eurasian Holarctic [past climate] and Australia)
ecospat.plot.contrib(contrib=pca.env.pastVau$co, eigen=pca.env.pastVau$eig) # plot variables contribution

# predict the scores on the axes
scores.globclim.pastVau <- pca.env.pastVau$li # PCA scores for the whole study area 
scores.sp.past.pastVau <- suprow(pca.env.pastVau, past.occ.original[which(past.occ.original[,"species_occ"]==1),bioclimsub])$li # PCA scores for occurrences in Eurasia (past climate)
head(scores.sp.past.pastVau)
scores.sp.present.pastVau <- suprow(pca.env.pastVau, Australia.occ[which(Australia.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for occurrences in Australia
scores.clim.past.pastVau <- suprow(pca.env.pastVau, past.occ.original[,bioclimsub])$li # PCA scores for all of the Eurasian Holarctic (past climate)
scores.clim.present.pastVau <- suprow(pca.env.pastVau, Australia.occ[,bioclimsub])$li # PCA scores for all of Australia

# Calculate the Occurrence Densities Grid (in environmental space, 100 x 100 cells)

# For the native range (past climate)
grid.clim.past.pastVau <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVau,
                                                glob1=scores.clim.past.pastVau, sp=scores.sp.past.pastVau, R=100, th.sp=0)

# For Australia
grid.clim.present.pastVau <- ecospat.grid.clim.dyn(glob=scores.globclim.pastVau,
                                                   glob1=scores.clim.present.pastVau, sp=scores.sp.present.pastVau, R=100, th.sp=0)  

# Calculate Niche Overlap (D value)
D.pastVau <- ecospat.niche.overlap (grid.clim.past.pastVau, grid.clim.present.pastVau, cor=TRUE)$D
D.pastVau

## Niche Equivalency Test
eq.test.pastVau <-ecospat.niche.equivalency.test(grid.clim.past.pastVau, grid.clim.present.pastVau,
                                                 rep=100, alternative = "greater") # setting alternative=greater, because D is higher than the null distribution
ecospat.plot.overlap.test(eq.test.pastVau, "D", "Equivalency")
eq.test.pastVau$p.D

## Niche Similarity Test
sim.test.pastVau <- ecospat.niche.similarity.test(grid.clim.past.pastVau, grid.clim.present.pastVau,
                                                  rep=100, alternative = "greater", rand.type = 2)  

ecospat.plot.overlap.test(sim.test.pastVau, "D", "Similarity")
sim.test.pastVau$p.D
# Result: niches are conserved between ranges (= more equivalent than expected by chance). 
# However, we cannot rule out that this is due to the availability of environments. 

# Delimiting niche categories and quantifying niche dynamics in analogue climates (expansion/stability/unfilling)
w.pastVau <- ecospat.niche.dyn.index(grid.clim.past.pastVau, grid.clim.present.pastVau, intersection=cutoff)$dynamic.index.w
w.pastVau 

# Visualizing niche categories, niche dynamics and climate analogy between ranges
myplot(grid.clim.past.pastVau, grid.clim.present.pastVau, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/blue
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
# Plot Variables Contribution
ecospat.plot.contrib(contrib=pca.env.pastVau$co, eigen=pca.env.pastVau$eig)

#### ================ PRESENT NATIVE RANGE CLIMATE AND OCCURRENCES TO AUSTRALIA =================== 

### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both the Eurasian Holarctic and Australia
pca.env.presentVau <- dudi.pca(rbind(present.occ, Australia.occ)[,bioclimsub], scannf=FALSE, nf=2)  #  this calibrates the PCA in the entire study area (Eurasia [present cliamte] and Australia)
ecospat.plot.contrib(contrib=pca.env.presentVau$co, eigen=pca.env.presentVau$eig) # plot variables contribution

# predict the scores on the axes
scores.globclim.presentVau <- pca.env.presentVau$li # PCA scores for the whole study area 
scores.sp.past.presentVau <- suprow(pca.env.presentVau, present.occ[which(present.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for the native occurrences (present)
head(scores.sp.past.presentVau)
scores.sp.present.presentVau <- suprow(pca.env.presentVau, Australia.occ[which(Australia.occ[,"species_occ"]==1),bioclimsub])$li # PCA scores for occurrences in Australia
scores.clim.past.presentVau <- suprow(pca.env.presentVau, present.occ[,bioclimsub])$li # PCA scores for the whole native study area (present climate) 
scores.clim.present.presentVau <- suprow(pca.env.presentVau, Australia.occ[,bioclimsub])$li # PCA scores for all of Australia

# Calculate the Occurrence Densities Grid (in environmental space, 100 x 100 cells)

# For the native range (present)
grid.clim.past.presentVau <- ecospat.grid.clim.dyn(glob=scores.globclim.presentVau,
                                                   glob1=scores.clim.past.presentVau, sp=scores.sp.past.presentVau, R=100, th.sp=0)

# For Australia
grid.clim.present.presentVau <- ecospat.grid.clim.dyn(glob=scores.globclim.presentVau,
                                                      glob1=scores.clim.present.presentVau, sp=scores.sp.present.presentVau, R=100, th.sp=0)  

# Calculate Niche Overlap (D value)
D.presentVau <- ecospat.niche.overlap (grid.clim.past.presentVau, grid.clim.present.presentVau, cor=TRUE)$D
D.presentVau # lower niche overlap when taking into account the native range expansion

# Perform the Niche Equivalency Test
eq.test.presentVau <-ecospat.niche.equivalency.test(grid.clim.past.presentVau, grid.clim.present.presentVau,
                                                    rep=100, alternative = "greater") 

ecospat.plot.overlap.test(eq.test.presentVau, "D", "Equivalency")
eq.test.presentVau$p.D

# Niche Similarity Test
sim.test.presentVau <- ecospat.niche.similarity.test(grid.clim.past.presentVau, grid.clim.present.presentVau,
                                                     rep=100, alternative = "greater", rand.type = 2)  

ecospat.plot.overlap.test(sim.test.presentVau, "D", "Similarity")
sim.test.presentVau$p.D

# Delimiting niche categories and quantifying niche dynamics in analogue climates (expansion/stability/unfilling)
w.presentVau <- ecospat.niche.dyn.index(grid.clim.past.presentVau, grid.clim.present.presentVau, intersection=cutoff)$dynamic.index.w
w.presentVau 

# Visualizing niche categories, niche dynamics and climate analogy between ranges
myplot(grid.clim.past.presentVau, grid.clim.present.presentVau, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "Niche Overlap", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/blue
       colZ1 = "#e29f60", colZ2 = "#52c6f5")

# Plot Variables Contribution
ecospat.plot.contrib(contrib=pca.env.presentVau$co, eigen=pca.env.presentVau$eig)


#### ================ FIGURE 3 ================

## Panel a-c
extAU <- extent(110,155,-45,-10) # spatial extent of Australia

#pdf("Figure3_top.pdf", width=8, height=2.8)
par(mfrow=c(1,3), mar=c(4, 4, 2, 1))

myplot(grid.clim.past.pastVau, grid.clim.present.pastVau, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/blue
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
myplot(grid.clim.past.presentVau, grid.clim.present.presentVau, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/blue
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
# Run LustenhouwerParker_MAxEnt.R first for the third panel, or skip these lines:
plot(crop(AUmap.difference.minPres, extAU), 
     col=c("grey30","pink","orange","grey90", "green",
           "red","grey50","yellow","purple","grey70"), main="", legend=F)
  points(y ~ x, Australia.records,cex=.4, pch=16) # black outline for points
  lines(Australia, lwd=1.5)
# end of third panel
dev.off()

## Panel d-f
USstates <- raster::getData('GADM', country='USA', level=1) 

#pdf("Figure3_bottom.pdf", width=8, height=2.8)
par(mfrow=c(1,3), mar=c(4, 4, 2, 1))

myplot(grid.clim.past.pastVca, grid.clim.present.pastVca, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/green        
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
myplot(grid.clim.past.presentVca, grid.clim.present.presentVca, quant=cutoff, 
       interest=1, # interest=1 plots native, interest=2 plots invasive density
       title= "", name.axis1="PC1", name.axis2="PC2",
       colz1 = "#eab98b", colz2 = "#9adef9", colinter = "#5f1988", # purple/orange/green
       colZ1 = "#e29f60", colZ2 = "#52c6f5")
# Run LustenhouwerParker_MAxEnt.R first for the third panel, or skip these lines:
plot(CAmap.difference.minPres, col=c("grey90","grey50","grey70"), 
     legend=F, main="")
  points(y ~ x, calflora.records, pch=16)
  lines(USstates)
# end of third panel
dev.off()

#### ================ TABLE 2 ===================

ecospatresults <- data.frame("niche1" = c("past", "past", "present", "past", "present"),
                             "niche2" = c("present", "Australia", "Australia", "California", "California"),
                             "D" = rep(NA,5),
                             "expansion" = rep(NA,5),
                             "stability" = rep(NA,5),
                             "unfilling" = rep(NA,5),
                             "H1.equivalency" = c("lower", "greater", "greater", "lower", "lower"),
                             "P.equivalency" = rep(NA,5),
                             "H1.similarity" = rep("greater",5),
                             "P.similarity" = rep(NA,5))
ecospatresults

ecospatresults$D <- c(D.pastVpresent, D.pastVau, D.presentVau, D.pastVca, D.presentVca)
ecospatresults[1,4:6] <- w.pastVpresent
ecospatresults[2,4:6] <- w.pastVau
ecospatresults[3,4:6] <- w.presentVau
ecospatresults[4,4:6] <- w.pastVca
ecospatresults[5,4:6] <- w.presentVca
  ecospatresults[,3:6] <-  round(ecospatresults[,3:6], 2)

ecospatresults[1,c(8,10)] <- c(eq.test.pastVpresent$p.D, sim.test.pastVpresent$p.D)
ecospatresults[2,c(8,10)] <- c(eq.test.pastVau$p.D, sim.test.pastVau$p.D)
ecospatresults[3,c(8,10)] <- c(eq.test.presentVau$p.D, sim.test.presentVau$p.D)
ecospatresults[4,c(8,10)] <- c(eq.test.pastVca$p.D, sim.test.pastVca$p.D)
ecospatresults[5,c(8,10)] <- c(eq.test.presentVca$p.D, sim.test.presentVca$p.D)
  ecospatresults[,c(8,10)] <- round(ecospatresults[,c(8,10)], 2)

ecospatresults # Table 2

