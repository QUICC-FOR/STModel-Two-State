#!/usr/bin/env Rscript

# depends:
#    dat/climateGrid_scaled.rds
#    dat/stm_calib/*
#    res/maps/*

# produces:
#    res/eval/ras/*

load("dat/map_projections.rdata")
source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')

climGrid = readRDS(file.path('dat', 'climateGrid_scaled.rds'))

subset_range = function(lon, lat, lonRange, latRange, spName = NA)
{
	if(spName == "183302-PIC-MAR") latLim[1] = latLim[1] + 1.9
	(lon > lonRange[1] &lon < lonRange[2] & lat > latRange[1] & lat < latRange[2])
}

suppressWarnings(dir.create('res/eval/ras', recursive = TRUE))



for(spName in speciesList)
{
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

	saveRDS(gr_env1, file.path('res', 'eval', 'ras', paste0(spName, '_env1.rds'))) 
	saveRDS(gr_env2, file.path('res', 'eval', 'ras', paste0(spName, '_env1.rds'))) 


	# sdm raster
	spGrid = readRDS(file.path('res','maps',paste0(spName,'_maps.rds')))
	spSubset = subset_range(spGrid$lon, spGrid$lat, lonLim, latLim, spName)
	sdmPres = make_raster(spGrid[spSubset,'sdm.pres'], 
			spGrid[spSubset, c('lon', 'lat')], P4S.latlon, stmMapProjection)
	stm = make_raster(spGrid[spSubset,'stm'], 
			spGrid[spSubset, c('lon', 'lat')], P4S.latlon, stmMapProjection)

	saveRDS(sdmPres, file.path('res', 'eval', 'ras', paste0(spName, '_sdmPres.rds'))) 
	saveRDS(stm, file.path('res', 'eval', 'ras', paste0(spName, '_stmProb.rds'))) 
}