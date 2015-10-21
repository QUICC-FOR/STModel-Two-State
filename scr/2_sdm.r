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
##       img/sdms.png

library(randomForest)
library(pROC)

# some constants
overwrite = FALSE
sdm.threshold = 0.25

dir.create(file.path('res', 'sdm'), recursive=TRUE)
dir.create(file.path('img'), recursive=TRUE)


arg = commandArgs(trailingOnly = TRUE)
if('--overwrite' %in% arg | '-o' %in% arg) overwrite = TRUE

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/plotClimate_scaled.rds")
climGrid = readRDS('dat/climateGrid_scaled.rds')
source("scr/stm_functions.r")

# subset the data for the SDM
# select only a single row for each plot (to avoid too much spatial duplication)
# note that all species get the same subset (will be applied during the species loop)
cat('Starting\nSelecting SDM rows\n')
rows = sapply(unique(climDat$plot_id), function(i) {
	candidates = which(climDat$plot_id == i)
	if(length(candidates) == 1) candidates else sample(candidates, 1)
	})


get_sdm_dat = function(sp, sel, climate)
{
	# select the data for a particular species to be used with the sdm
	# sp: the species code
	# baseDir: the directory containing all the species' data
	# sel: the rows in the dataframe to select for the SDM
	# climate: the climate data to merge with the SDM data
	# returns a list with two dataframes
	# selected is the dataset for use in the SDM
	# unselected is the rest of the data
	
	spPath = file.path('dat', 'presence', paste0(sp, '_presence.rds'))
	spData = readRDS(spPath)
	colnames(spData)[3] = 'presence'
	sdmDat = merge(spData, climate)

	# apply the subset, and drop the plot ID and year measured columns also
	sdmDat.subset = sdmDat[sel,-c(1,2)]
	sdmDat.unsubset = sdmDat[-sel,-c(1,2)]

	list(selected = sdmDat.subset, unselected = sdmDat.unsubset)
}


get_auc = function(mod, newdata, presence.name = 'presence')
{
	## computes auc for the rf model
	# mod: the rf model to use
	# newdata: a dataframe containing independent data for the evaluation
	# presence.name: the name or index of the column containing the presence data
	# see how well the rf fits
	unsel.preds = predict(mod, newdata = newdata, type='prob')[,2]
	rf.roc = roc(response=newdata[,presence.name], predictor=unsel.preds)
	return(rf.roc$auc)
}


sdm_pres = function(sdm, newdata, threshold)
{
	# return a presence absence vector using the predictions of an sdm and a threshold
	predictions = predict(sdm, newdata = newdata, type='prob')[,2]
	as.integer(predictions >= threshold)

}




#######
# MAIN LOOP ACROSS SPECIES
#######

rocResults = list()
sdmModels = list()
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
	} else {
		sdmDat = get_sdm_dat(spName, rows, climDat)
		
		## run RF
		cat("  Fitting random forest\n")
		rf.mod = randomForest(as.factor(presence) ~ . , data = sdmDat$selected, ntree = 500)
		saveRDS(rf.mod, rfModFilename)
		sdmModels[[spName]] = rf.mod
		
		# project SDM into geographic space
		cat("  Projecting SDM geographically\n")
		sdmProjection = data.frame(lon=climGrid$lon, lat=climGrid$lat, 
				sdm=predict(rf.mod, newdata=climGrid, type='prob')[,2])
		sdmProjection$sdm.pres = sdm_pres(rf.mod, climGrid, sdm.threshold)
## 		sdmProjection$stm.mask = with(sdmProjection, stm_mask(cbind(lon, lat), 
## 				presDat[,spName], presDat[,c('lon', 'lat')], stmMaskTolerance))
		saveRDS(sdmProjection, projFilename)
		projections[[spName]] = sdmProjection

		
		cat("  Computing ROC\n")
		rocResults[[spName]] = data.frame(species=spName, roc=get_auc(rf.mod, sdmDat$unselected))
	}
}


# write ROC
rocTable = do.call(rbind, rocResults)
write.csv(rocTable, file=file.path('res', 'sdmROC.csv'))
	
# plot the SDMs
coord.cols=c('lon', 'lat')
fname = file.path('img', 'sdms.png')
set_up_stm_figure(length(speciesList) + 1, fname, mar=c(0,0,1,0), oma=c(0.5,0.5,2.5,0.5))
for(spName in speciesList)
	plot_sdm(projections[[spName]]$sdm, climGrid[,coord.cols], sdmColors, zlim=c(0,1), main=spName)

# plot a legend
plot.new()
par(mar=c(15,0,0,0))
plot_sdm(projections[[spName]]$sdm, climGrid[,coord.cols], sdmColors, zlim=c(0,1), legend=TRUE, 
		legend.only=TRUE, plot.ocean=FALSE, legend.width=6, legend.shrink=1, horizontal=TRUE)
clean_up_figure(fname)

