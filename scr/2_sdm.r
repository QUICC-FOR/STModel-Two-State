#!/usr/bin/Rscript
library(randomForest)
library(pROC)

# some constants
overwrite = FALSE


arg = commandArgs(trailingOnly = TRUE)
if('--overwrite' %in% arg | '-o' %in% arg) overwrite = TRUE

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/plotClimate_scaled.rds")

# subset the data for the SDM
# select only a single row for each plot (to avoid too much spatial duplication)
# note that all species get the same subset (will be applied during the species loop)
cat('Starting\nSelecting SDM rows\n')
rows = sapply(unique(climDat$plot_id), function(i) {
	candidates = which(climDat$plot_id == i)
	if(length(candidates) == 1) candidates else sample(candidates, 1)
	})


get_sdm_dat = function(sp, baseDir, sel, climate)
{
	# select the data for a particular species to be used with the sdm
	# sp: the species code
	# baseDir: the directory containing all the species' data
	# sel: the rows in the dataframe to select for the SDM
	# climate: the climate data to merge with the SDM data
	# returns a list with two dataframes
	# selected is the dataset for use in the SDM
	# unselected is the rest of the data
	
	spPath = file.path(baseDir, 'dat', paste(sp, 'presence.rds', sep='_'))
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


#######
# MAIN LOOP ACROSS SPECIES
#######

for(spName in speciesList)
{
	cat(paste("Starting species", spName, '\n'))

	baseDir = file.path('species', spName)
	rfModFilename = file.path(baseDir, 'res', paste(spName, 'sdm.rds', sep='_'))

	if(file.exists(rfModFilename) & !overwrite)
	{
		warning(paste("Output file", rfModFilename, 
				"already exists; skipping. Use --overwrite or -o to redo all species"))
	} else {
		sdmDat = get_sdm_dat(spName, baseDir, rows, climDat)
		
		## run RF
		cat("  Fitting random forest\n")
		rf.mod = randomForest(as.factor(presence) ~ . , data = sdmDat$selected, ntree = 500)
		saveRDS(rf.mod, rfModFilename)
		
		cat("  Computing ROC\n")
		cat("Area under the ROC curve: ", get_auc(rf.mod, sdmDat$unselected), "\n", 
			file=file.path(baseDir, 'res', 'sdm_roc.txt'))
		}
}