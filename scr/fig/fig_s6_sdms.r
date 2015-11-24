#!/usr/bin/env Rscript

## produce a plot of the SDMs
## dependencies:
##       dat/speciesList.rds  (script 1)
##       dat/map_projections.rds  (script 1)
##       dat/ne_50m_ocean
##       dat/little_rangemaps/*
##       res/sdm/*_sdm_projection.rds

## creates:
##       img/figs/fig_s6.pdf

library(randomForest)
library(rgdal)
library(sp)

speciesList = readRDS('dat/speciesList.rds')

draw.plots = TRUE
arg = commandArgs(trailingOnly = TRUE)
if('--plots' %in% arg) draw.plots = TRUE
if('--noplots' %in% arg) draw.plots = FALSE

source("scr/stm_functions.r")
climDat = readRDS("dat/clim/plotClimate_scaled.rds")

speciesInfo = read.csv('dat/speciesInfo.csv')
if(draw.plots)
{
	plotLocs = readRDS("dat/raw/plotLocations.rds")
	climDat = merge(climDat, plotLocs)
}

cex.ti = 1.5
figure.width = 16
figure.height = 24
fname = file.path('img', 'figs', 'fig_s6.png')
png(file=fname, width=figure.width, height=figure.height, pointsize=9, res=300, units='cm')
## pdf(file=fname, width=figure.width/2.54, height=figure.height/2.54, pointsize=9)
par(mfrow=c(6, 4), mar=c(0,0,2,0), oma=c(1,1,0,2), tcl=-0.2)

cols = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"), 
		interpolate='spline', bias=1, space="rgb")(200)

for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	sdmFilename = file.path('res', 'sdm', paste0(spName, "_sdm_projection.rds"))
	sdm.projection = readRDS(sdmFilename)
	plot.sdm(sdm.projection, cols, axes=FALSE, box=TRUE, main=plLab)
## 	mtext(plLab, side=3)
	if(draw.plots)
	{
		spPath = file.path('dat', 'presence', paste0(spName, '_presence.rds'))
		spData = readRDS(spPath)
		colnames(spData)[3] = 'presence'
		sdmDat = merge(spData, climDat)
		sdmDat = sdmDat[sdmDat$presence == 1,]
		sdmDat = sdmDat[sample(nrow(sdmDat), as.integer(0.1*nrow(sdmDat))),]
		with(sdmDat, points(x, y, pch=16, cex=.2, col='#000000cc'))
	}
}


# plot a legend
plot.new()
par(mar=c(15,0,0,0))
plot.sdm(sdm.projection, cols, zlim=c(0,1), legend=TRUE, legend.only=TRUE, 
		plot.ocean=FALSE, legend.width=6, legend.shrink=1, horizontal=TRUE, smallplot=c(0.1, 0.9, 0.2, 0.35), legend.cex=1)

plot.new()
par(xpd=NA)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(temperate.ypos, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(transitional.ypos, "ndc", "user"), "Transitional", 
		srt=270, cex=cex.ti)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(boreal.ypos, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)

dev.off()