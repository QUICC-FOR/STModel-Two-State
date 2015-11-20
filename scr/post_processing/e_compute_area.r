#!/usr/bin/env Rscript

## depends:
##    res/posteriorGrid/*

## makes:
##    res/maps/*

speciesList = readRDS('dat/speciesList.rds')
models = c('0', 'i0', 'g', 'ig')

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

for(spName in speciesList)
{
	# get the calibration range
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
	latLim = range(calibDat$lat)
	lonLim = range(calibDat$lon)

	for(mod in models)
	{
		spGrid = readRDS(file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
		postGrid = readRDS(file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
		grPres = postGrid$stmPres
		lon = spGrid$lon
		lat = spGrid$lat

		areas = matrix(NA, nrow=nrow(grPres), ncol=3)
		for(i in 1:nrow(grPres))
		{
			pres = grPres[i,]
			rde = (1 * (pres & spGrid$sdm.pres)) + (2 * (pres & !spGrid$sdm.pres)) + 
						(3 * (!pres & spGrid$sdm.pres))
			rde[rde == 0] = NA	
			rde = rde - 1
			
		# restrict to calibration range
			rde[lon < lonLim[1] | lon > lonLim[2] | lat < latLim[1] | lat > latLim[2]] = NA
			if(spName == "183302-PIC-MAR")
				rde[lat < (latLim[1]+1.9)] = NA

			rdeRas = make_raster(rde, spGrid[,1:2], P4S.latlon, stmMapProjection)
			areas[i,] = freq(rdeRas)[1:3,2] * prod(res(rdeRas)/1000)/1000
		}
		areas = as.data.frame(areas)
		colnames(areas) = c('present', 'expand', 'contract')
		saveRDS(areas, file.path('res','areas',paste0(spName, '_', mod, '_areas.rds')))
	}	
}