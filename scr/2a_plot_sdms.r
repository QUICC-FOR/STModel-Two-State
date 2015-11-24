#!/usr/bin/env Rscript

## produce a plot of the SDMs
## dependencies:
##       dat/speciesList.rds  (script 1)
##       dat/map_projections.rds  (script 1)
##       dat/ne_50m_ocean
##       dat/little_rangemaps/*
##       res/sdm/*_sdm_projection.rds

## creates:
##       img/sdms.png

library(randomForest)
library(rgdal)
library(sp)

speciesList = readRDS('dat/speciesList.rds')
dir.create(file.path('img'), recursive=TRUE)

draw.plots = FALSE
arg = commandArgs(trailingOnly = TRUE)
if('--plots' %in% arg) plots = TRUE

source("scr/stm_functions.r")
climDat = readRDS("dat/clim/plotClimate_scaled.rds")
plotLocs = readRDS("dat/raw/plotLocations.rds")
climDat = merge(climDat, plotLocs)

# plot the SDMs

fname = file.path('img', 'sdms.png')
set_up_stm_figure(length(speciesList) + 1, fname, mar=c(0,0,1,0), oma=c(0.5,0.5,2.5,0.5))
for(spName in speciesList)
{
	sdmFilename = file.path('res', 'sdm', paste0(spName, "_sdm_projection.rds"))
	sdm.projection = readRDS(sdmFilename)
	plot.sdm(sdm.projection, sdmColors, main=spName)
	if(file.exists(paste0('dat/little_rangemaps/', spName)))
	{
		range.map = readOGR(dsn=paste0('dat/little_rangemaps/', spName), layer=spName)
		proj4string(range.map) = CRS("+proj=longlat +datum=NAD27")
		range.map = spTransform(range.map, mapProj$projected)
		plot(range.map, border='red', col="#FFFFFF00", add=TRUE)

	}
	if(draw.plots)
	{
		spPath = file.path('dat', 'presence', paste0(spName, '_presence.rds'))
		spData = readRDS(spPath)
		colnames(spData)[3] = 'presence'
		sdmDat = merge(spData, climDat)
		with(sdmDat[sdmDat$presence == 1,], points(x, y, pch=3, cex=.3, col='#ffaaaa22'))
	}
}


# plot a legend
## plot.new()
## par(mar=c(15,0,0,0))
## plot.sdm(projections[[spName]], sdmColors, zlim=c(0,1), legend=TRUE, legend.only=TRUE, 
## 		plot.ocean=FALSE, legend.width=6, legend.shrink=1, horizontal=TRUE)
clean_up_figure(fname)

# plot all the climate variables
library(raster)
climGrid = readRDS('dat/clim/climateGrid_unscaled.rds')
coordinates(climGrid) = c('x', 'y')
gridded(climGrid) = TRUE
climGrid = brick(climGrid)

dpi = 600
figure.width = 12
figure.height = 12
filename = file.path('img', 'sdm_climvars.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)

par(mfrow=c(2,3))
for(cv in names(climGrid)) {
	plot(climGrid[[cv]], legend=FALSE, main=cv)
	plot(ocean, col="white", add=TRUE)
	if(draw.plots) (with(sdmDat, points(x,y, pch=3, cex=0.2, col='#00000011')))
}
dev.off()