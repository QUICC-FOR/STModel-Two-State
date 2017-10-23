#!/usr/bin/env Rscript


library("raster")
library("INLA")
library('rgdal')
source("scr/inla/inla_functions.r")

spAll <- readRDS("dat/speciesList.rds")
spList <- commandArgs()
spList <- spList[spList %in% spAll]
if(length(spList) == 0) spList <- spAll

climGrid <- readRDS("dat/clim/climateGrid_scaled.rds")
basedat <- readRDS('dat/basedata.rds')

for(spName in spList)
{
	cat(paste0(spName, '\n'))
	# get models
	resDir <- file.path('res', 'inla', spName)
	mods <- readRDS(file.path(resDir, 'inla_models_ext.rds'))


	## RANDOM FIELDS
	imgDir <- 'img/inla/gp_plots/'
	dir.create(imgDir, showWarnings=FALSE, recursive=TRUE)
	fname <- paste0(spName, '_gfields.pdf')
	pdf(file.path(imgDir, fname), w=12, h=12)
	par(mfrow=c(1,2), mar=c(0.5,0.5,3,0.5), oma=c(0,0,0,3))
	# extinction model
	extGP <- inla_random_field(mods$mesh, mods$models[[2]], range(climGrid$x), range(climGrid$y))
	inla_rf_plots(extGP, basedata=basedat, main="Extinction Gaussian Field")

	dead <- dev.off()
}








