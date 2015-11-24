#!/usr/bin/env Rscript

## depends:
##    res/posteriorGrid/*

## makes:
##    res/areas/*

library(foreach)
library(doParallel)
numCores = detectCores() - 1
registerDoParallel(cores=numCores)

speciesList = readRDS('dat/speciesList.rds')
mod = '0'

arg = commandArgs(TRUE)
if(length(arg) > 0)
{
	sps = which(speciesList %in% arg)
	if(length(sps) == 0) 
	{
		warning("No species specified on command line; falling back to default")
	} else 
		speciesList = speciesList[which(speciesList %in% arg)]
}

if(length(speciesList) == 0) stop("Error: no species specified")
cat("Will process posteriors for:\n")
for(sp in speciesList) cat("  ", sp, "\n")

source('scr/stm_functions.r')
load("dat/map_projections.rdata")
suppressWarnings(dir.create(file.path('res', 'areas'), recursive=TRUE))


cat('starting species loop\n')
for(spName in speciesList)
{
	cat("  ", spName, "\n")
	cat('    reading calibration data\n')
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
	latLim = range(calibDat[,'lat'])
	lonLim = range(calibDat[,'lon'])

	cat('    reading grids\n')	
	spGrid = readRDS(file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
	postGrid = readRDS(file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
	grPres = postGrid$stmPres
	lon = spGrid$lon
	lat = spGrid$lat

	cat('    computing across posterior\n')	
	area_sp = foreach(pres = iter(grPres, by='row'), .combine=rbind, 
	.packages=c('raster', 'rgdal'), .final=function(x) {
			colnames(x) = c('present', 'expand', 'contract')
			as.data.frame(x)
		}) %dopar% {
			rde = (1 * (pres & spGrid$sdm.pres)) + (2 * (pres & !spGrid$sdm.pres)) + 
						(3 * (!pres & spGrid$sdm.pres))
			rde[rde == 0] = NA	
			rde = rde - 1
			
			# restrict to calibration range
			rde[lon < lonLim[[spName]][1] | lon > lonLim[[spName]][2] | lat < latLim[[spName]][1] | lat > latLim[[spName]][2]] = NA
			if(spName == "183302-PIC-MAR")
				rde[lat < (latLim[[spName]][1]+1.9)] = NA

			rdeRas = make_raster(rde, spGrid[,1:2], P4S.latlon, stmMapProjection)
			freq(rdeRas)[1:3,2] * prod(res(rdeRas)/1000)/1000
	}
	cat('    saving\n')	
	saveRDS(area_sp, file.path('res','areas',paste0(spName, '_', mod, '_areas.rds')))
}
cat('done\n')
