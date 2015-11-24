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
models = c('0', 'g')

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


# get calibration data
cat('reading calibration data\n')
calibDat = sapply(speciesList, simplify = FALSE, USE.NAMES = TRUE, FUN=function(spName) 
		readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds'))))	
	
# compute lat/lon limits
latLim = lapply(calibDat, FUN=function(x) range(x[,'lat']))
lonLim = lapply(calibDat, function(x) range(x[,'lon']))


cat('starting species loop\n')
areas = foreach(spName = speciesList, .final = function(x) {names(x) = speciesList; x}) %dopar%
{
	foreach(mod = models, .final = function(x) {names(x) = models; x}) %do%
	{
		spGrid = readRDS(file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
		postGrid = readRDS(file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
		grPres = postGrid$stmPres
		lon = spGrid$lon
		lat = spGrid$lat

		area_sp = matrix(NA, nrow=nrow(grPres), ncol=3)
		for(i in 1:nrow(grPres))
		{
			pres = grPres[i,]
			rde = (1 * (pres & spGrid$sdm.pres)) + (2 * (pres & !spGrid$sdm.pres)) + 
						(3 * (!pres & spGrid$sdm.pres))
			rde[rde == 0] = NA	
			rde = rde - 1
			
		# restrict to calibration range
			rde[lon < lonLim[[spName]][1] | lon > lonLim[[spName]][2] | lat < latLim[[spName]][1] | lat > latLim[[spName]][2]] = NA
			if(spName == "183302-PIC-MAR")
				rde[lat < (latLim[[spName]][1]+1.9)] = NA

			rdeRas = make_raster(rde, spGrid[,1:2], P4S.latlon, stmMapProjection)
			area_sp[i,] = freq(rdeRas)[1:3,2] * prod(res(rdeRas)/1000)/1000
		}
		area_sp = as.data.frame(area_sp)
		colnames(area_sp) = c('present', 'expand', 'contract')
		area_sp
	}	
}

cat('saving files\n')
for(spName in speciesList) {
	for(mod in models) {
		saveRDS(areas[[spName]][[mod]], 
				file.path('res','areas',paste0(spName, '_', mod, '_areas.rds')))
	}
}
cat('done\n')
