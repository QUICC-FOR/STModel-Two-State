#!/usr/bin/env Rscript

# depends:
#    res/posterior/*
#    dat/climateGrid_scaled.rds
#    dat/stm_calib/*
#    res/maps/*

# produces:
#    res/eval/SPNAME_MOD__rdeAreas.rds

library(foreach)
library(doParallel)
library(coda)
library(raster)
source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')

suppressWarnings(dir.create('res/eval', recursive = TRUE))

load("dat/map_projections.rdata")

clArgs <- commandArgs(trailingOnly = TRUE)

## these will be read in from the command line
spName = clArgs[1]
modName = if(length(clArgs) > 1) clArgs[2] else '0'
numCores = if(length(clArgs) > 2) clArgs[3] else detectCores()
sampleSize = if(length(clArgs) > 3) clArgs[4] else NA

registerDoParallel(cores=numCores)


subset_range = function(lon, lat, lonRange, latRange, spName = NA)
{
	if(spName == "183302-PIC-MAR") latLim[1] = latLim[1] + 1.9
	(lon > lonRange[1] &lon < lonRange[2] & lat > latRange[1] & lat < latRange[2])
}


climGrid = readRDS(file.path('dat', 'climateGrid_scaled.rds'))


# get the calibration range
calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
latLim = range(calibDat$lat)
lonLim = range(calibDat$lon)

#create 2 projected rasters; env1 and env2
grSubset = subset_range(climGrid$lon, climGrid$lat, lonLim, latLim, spName)
gr_env1 = make_raster(climGrid[grSubset,'annual_mean_temp'], 
		climGrid[grSubset, c('lon', 'lat')], P4S.latlon, stmMapProjection)
gr_env2 = make_raster(climGrid[grSubset,'tot_annual_pp'], 
		climGrid[grSubset, c('lon', 'lat')], P4S.latlon, stmMapProjection)


# sdm raster
spGrid = readRDS(file.path('res','maps',paste0(spName,'_maps.rds')))
spSubset = subset_range(spGrid$lon, spGrid$lat, lonLim, latLim, spName)
sdmPres = make_raster(spGrid[spSubset,'sdm.pres'], 
		spGrid[spSubset, c('lon', 'lat')], P4S.latlon, stmMapProjection)

rdeCats = list('absent' = 0, 'present' = 1, 'expand' = 2, 'contract' = 3)

# grab the posterior samples
posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
posterior = readRDS(posteriorFname)[[modName]]

if(is.null(sampleSize) || is.na(sampleSize)) {
	samples = posterior
} else
	samples = posterior[sample(nrow(posterior), sampleSize),]


# foreach posterior replicate
rdePosterior = foreach(curPars = iter(samples, by='row'), .combine=rbind, 
				.packages = 'raster', .final=function(x) {
	colnames(x) = names(rdeCats)[2:4]
	x}) %dopar%
{
	# 1. compute C & E (as rasters! -- this saves having to reproject at each rep)
	if(length(curPars) == 2)
	{
		cPars = curPars[1]
		ePars = curPars[2]
	} else {
		cPars = curPars[1:5]
		ePars = curPars[6:10]
	}
	C = predict.stm_point(cPars, gr_env1, gr_env2)
	E = predict.stm_point(ePars, gr_env1, gr_env2)

	# 2. compute STM presence
	stmPres = (C > E)

	# 3. compute RDE
	rde = (rdeCats[['present']] * (stmPres & sdmPres)) + 
		(rdeCats[['expand']] * (stmPres & !sdmPres)) + 
		(rdeCats[['contract']] * (!stmPres & sdmPres))
	# change absences to NA and subtract 1 from all cats
	rde[rde==0] = NA
	rde = rde - 1
	
	# 5. compute zonal area
	cellArea_km2 = prod(res(rde)/1000)
	freqs = freq(rde)
	(freqs[,2] * cellArea_km2)[1:3]
}

outFile = file.path('res','eval',paste0(spName, '_', modName, '_rdeAreas.rds'))
saveRDS(rdePosterior, outFile)