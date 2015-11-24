#!/usr/bin/env Rscript

## perform SDM on all species
## dependencies:
##       dat/speciesList.rds  (script 1)
##       dat/plotClimate_scaled.rds  (script 1)
##       dat/climGrid_scaled.rds  (script 1)
##       dat/presence/*

## creates:
##       res/sdm/*
##       res/sdmROC.csv

library(randomForest)
library(pROC)

# some constants
overwrite = FALSE

dir.create(file.path('res', 'sdm'), recursive=TRUE)


arg = commandArgs(trailingOnly = TRUE)
if('--overwrite' %in% arg | '-o' %in% arg) overwrite = TRUE

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/clim/plotClimate_scaled.rds")
climGrid = readRDS('dat/clim/climateGrid_scaled.rds')



#######
# MAIN LOOP ACROSS SPECIES
#######

projections = list()
for(spName in speciesList)
{
	cat(paste("Starting species", spName, '\n'))
	rfModFilename = file.path('res', 'sdm', paste0(spName, "_rf_sdm.rds"))
	projFilename = file.path('res', 'sdm', paste0(spName, "_sdm_projection.rds"))

	if(file.exists(rfModFilename) & !overwrite)
	{
		warning(paste("Output file", rfModFilename, 
				"already exists; skipping. Use --overwrite or -o to redo all species"))
		rf.mod = readRDS(rfModFilename)
	} else {
		spPath = file.path('dat', 'presence', paste0(spName, '_presence.rds'))
		spData = readRDS(spPath)
		colnames(spData)[3] = 'presence'
		sdmDat = merge(spData, climDat)[,-c(1,2)]
		
		## run RF
		cat("  Fitting random forest\n")
		rf.mod = randomForest(as.factor(presence) ~ . , data = sdmDat, ntree = 500)
		saveRDS(rf.mod, rfModFilename)
	}
		
	# project SDM into geographic space
	cat("  Projecting SDM geographically\n")
	sdmProjection = data.frame(x=climGrid$x, y=climGrid$y, 
			sdm=predict(rf.mod, newdata=climGrid, type='prob')[,2])
	saveRDS(sdmProjection, projFilename)
	projections[[spName]] = sdmProjection
}



