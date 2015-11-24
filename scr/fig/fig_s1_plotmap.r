#!/usr/bin/env Rscript

## 
## source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')
## speciesInfo = read.csv('dat/speciesInfo.csv')
## suppressWarnings(dir.create(file.path('img', 'figs'), recursive=TRUE))
## 

library(sp)
library(raster)


# basedata
basedata = readRDS("dat/basedata.rds")

# grab the calib data from one species
calDat = readRDS(paste0('dat/transition/', speciesList[1], '_transitions.rds'))
plotLocs = readRDS('dat/raw/plotLocations.rds')
calDat = merge(calDat, plotLocs, by.x='plot', by.y='plot_id')

# reduce to unique plots only (since plots are duplicated) and make spatial
calDat = calDat[,c('x', 'y', 'plot')]
calDat = unique(calDat)
coordinates(calDat) = c('x', 'y')

# get the validation points
evalDat = readRDS(paste0('dat/eval/', speciesList[1], '_evalCells.rds'))
evalDat = evalDat[,c('x', 'y')]
evalDat = unique(evalDat)
coordinates(evalDat) = c('x','y')

water.col = "#7EB2D7"
ocean.border="#444444"
lake.border = "#007FD400"
river.color = "#007FD4"
plot.color = '#000000'
subsample.size = 0.1

dpi = 600
figure.width = 6.5
figure.height= 4
filename = file.path('img', 'figs', 'fig_s1.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)

par(mfrow=c(1,2), mar=c(0,0,0,0.5), oma=c(0.5,0.5,0.5,0))
plot(calDat, cex=0)
## plotRGB(basedata$hillshade, add=TRUE)
## plotRGB(basedata$nat.earth1, add=TRUE)
plotRGB(basedata$nat.earth2, add=TRUE)
points(calDat[sample(nrow(calDat), floor(subsample.size*nrow(calDat))),], cex=0.3, pch=3, col=plot.color)
## plotRGB(basedata$hypso, add=TRUE)
plot(basedata$rivers, add=TRUE, col=river.color)
## plot(basedata$ocean, add=TRUE, col=water.col, border=ocean.border)
plot(basedata$lakes, add=TRUE, col=water.col, border=lake.border)

text(grconvertX(0.05, "npc", "user"), grconvertY(0.95, "npc", "user"), "A")
## quartz()

plot(calDat, cex=0)
plotRGB(basedata$nat.earth2, add=TRUE)
points(evalDat[sample(length(evalDat), floor(subsample.size*length(evalDat))),], cex=0.3, pch=3, col=plot.color)
plot(basedata$rivers, add=TRUE, col=river.color)
plot(basedata$lakes, add=TRUE, col=water.col, border=lake.border)
text(grconvertX(0.05, "npc", "user"), grconvertY(0.95, "npc", "user"), "B")

dev.off()

