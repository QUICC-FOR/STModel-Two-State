#!/usr/bin/Rscript
library(randomForest)
library(pROC)
library(raster)
library(rgdal)

print("Setting up variables and reading files")
annealFraction = 0.5

speciesList = readRDS('dat/speciesList.rds')
transitionClimate = readRDS('dat/transitionClimate_scaled.rds')
load('dat/map_projections.rdata')

ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)

climDat = readRDS("dat/plotClimate_scaled.rds")


# subset the data for the SDM
# select only a single row for each plot (to avoid too much spatial duplication)
# note that all species get the same subset (will be applied during the species loop)
rows = sapply(unique(sdmDat$plot_id), function(i) {
	candidates = which(sdmDat$plot_id == i)
	if(length(candidates) == 1) candidates else sample(candidates, 1)
	})


for(spName in speciesList)
{
	print(paste("Starting species", spName))
	baseDir = file.path('species', spName)

	spPath = file.path(baseDir, 'dat', paste(spName, 'presence.rds', sep='_'))
	spData = readRDS(spPath)
	colnames(spData)[3] = 'presence'
	sdmDat = merge(spData, climDat)

	# apply the subset, and drop the plot ID and year measured columns also
	sdmDat.subset = sdmDat[rows,-c(1,2)]
	sdmDat.unsubset = sdmDat[-rows,-c(1,2)]

	print("  Fitting RF model")
	rf.mod = randomForest(as.factor(presence) ~ . , data = sdmDat.subset, ntree = 500)
	saveRDS(rf.mod, file.path(baseDir, 'res', paste(spName, 'sdm.rds', sep='_'))

	# see how well the rf fits
	print("  Computing ROC")
	ubsubset.predictions = predict(rf.mod, newdata = sdmDat.unsubset, type='prob')[,2]
	rf.roc = roc(response=sdmDat.unsubset$presence, predictor=ubsubset.predictions)
	cat("Area under the ROC curve: ", rf.roc$auc, "\n", 
			file=file.path(baseDir, 'res', 'sdm_roc.txt'))

	# project the SDM prob of presence to the transition observations
	print("  Projecting SDM to transitions")
	transitionData = readRDS(file.path(baseDir, 'dat', paste(spName, 'transitions.rds', sep='_')))
	stmData = merge(transitionData, transitionClimate)
	stmData$prevalence = predict(rf.mod, newdata = stmData, type='prob')[,2]
	
	# project SDM onto the climate grid
	print("  Projecting SDM geographically")
	climGrid = readRDS('dat/climateGrid_scaled.rds')
	sdmProjection = data.frame(lon=climGrid$lon, lat=climGrid$lat, 
			sdm=predict(rf.mod, newdata=climGrid, type='prob')[,2])
	saveRDS(sdmProjection, file.path(baseDir, 'res', paste(spName, 'sdm_grid_projection.rds', sep='_')))


	# subset the transition data for use in the annealing
	sel = sample(nrow(stmData), as.integer(annealFraction*nrow(stmData)))
	stmData.subset = stmData[sel,]
	stmData.unsubset = stmData[-sel,]
	saveRDS(stmData.subset, file.path(baseDir, 'dat', paste(spName, 'stm', 'calib.rds', sep='_')))
	saveRDS(stmData.unsubset, file.path(baseDir, 'dat', paste(spName, 'stm', 'valid.rds', sep='_')))
	
	# diagnostic plots
	# response curve
	print("  Producing SDM plots")
	pdf(file=file.path(baseDir, 'img', 'sdm_response.pdf'), width=10,height=7)
	par(mfrow=c(2,3), bty='n', oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))
	xx = matrix(0, nrow=1000, ncol=6)
	xx = as.data.frame(xx)
	colnames(xx) = colnames(sdmDat)[-c(1:3)]
	for(i in 1:ncol(xx))
	{
		xx[,i] = seq(-3, 3, length.out=1000)
		yy = predict(rf.mod, newdata=xx, type='prob')[,2]
		plot(xx[,i], yy, type='l', col='blue', xlab=colnames(xx)[i], ylab="prob of presence", ylim=c(0,1))
		xx[,i] = rep(0,1000)
	}
	mtext(spName, outer=T, line=0.7, cex=0.8)
	dev.off()
	
	# map
	coordinates(sdmProjection) = c('lon', 'lat')
	gridded(sdmProjection) = TRUE
	sdmProjection = raster(sdmProjection)
	load('dat/map_projections.rdata')
	proj4string(sdmProjection) = P4S.latlon
	sdmProjection = projectRaster(sdmProjection, crs=stmMapProjection)
	
	sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
			interpolate='spline', bias=1, space="rgb")(200)
	pdf(file=file.path(baseDir, 'img', 'sdm.pdf'), width=7, height=7)
	plot(sdmProjection, main=paste("SDM:", spName), xaxt='n', yaxt='n', col=sdmColors)
	plot(ocean, col="white", add=TRUE)
	dev.off()

}


