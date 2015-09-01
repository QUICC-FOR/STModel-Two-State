#!/usr/bin/Rscript
library(randomForest)
library(pROC)
library(raster)
library(rgdal)

# some constants
overwrite = FALSE
sdm.threshold = 0.25
annealFrac = 0.5
stmMaxInterval = 15
mask.tol = 2 # lat/lon tolerance for the data mask (how many degrees out of the range for projecting/calibrating
sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
		interpolate='spline', bias=1, space="rgb")(200)


arg = commandArgs(trailingOnly = TRUE)
if('--overwrite' %in% arg | '-o' %in% arg) overwrite = TRUE

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/plotClimate_scaled.rds")
transitionClimate = readRDS('dat/transitionClimate_scaled.rds')
climGrid = readRDS('dat/climateGrid_scaled.rds')
plotLocs = readRDS('dat/plotLocations.rds')
transitionClimate = merge(transitionClimate, plotLocs, by.x=1, by.y=1)
load('dat/map_projections.rdata')
ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)

# subset the data for the SDM
# select only a single row for each plot (to avoid too much spatial duplication)
# note that all species get the same subset (will be applied during the species loop)
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


make_raster = function(dat, coords, start.proj = NULL, dest.proj = NULL)
{
	# simple utility to make a projected raster out of a vector and lat/long coords
	require(raster)
	require(rgdal)
	if(!(ncol(coords) == 2 & nrow(coords) == length(dat)))
		stop("Coords must have 2 columns and a number of rows equal to the length of dat")
	ras = cbind(dat, coords)
	coordinates(ras) = c(2,3)
	gridded(ras) = TRUE
	ras = raster(ras)
	if(!is.null(start.proj))
	{
		proj4string(ras) = P4S.latlon
		if(!is.null(dest.proj))
			ras = projectRaster(ras, crs=dest.proj)
	}
	return(ras)
}




compute_mask = function(sdm, newdata, threshold, tol = c(10,10))
{
	# computes the species range from an sdm (a random forest)
	# returns data frame with 2 columns; the first a presence-absence vector for the sdm
	# and the second a mask vector showing the data limits within tol degrees of lon/lat
	# if tol is length 1, it will be applied in both dimension
	# if length 2, then the first will be longitude, second latitude

	if(length(tol) == 1) tol = rep(tol,2)

	# first compute the sdm presence
	predictions = predict(sdm, newdata = newdata, type='prob')[,2]
	sdm.pres = as.integer(predictions >= threshold)
	
	# next get the mask
	lon.r = with(newdata, range(lon[sdm.pres == 1])) + tol[1]*c(-1,1)
	lat.r = with(newdata, range(lat[sdm.pres == 1])) + tol[2]*c(-1,1)
	mask.ind = with(newdata, 
			which(lon >= lon.r[1] & lon <= lon.r[2] & lat >= lat.r[1] & lat <= lat.r[2]))
	sdm.mask = rep(0, nrow(newdata))
	sdm.mask[mask.ind] = 1
	
	data.frame(sdm.pres, sdm.mask)
}


prep_transitions = function(sp, baseDir, climate, sdm, thresh, tol)
{
	transitionData = readRDS(file.path(baseDir, 'dat', 
			paste(sp, 'transitions.rds', sep='_')))
	prevalence = predict(sdm, newdata = climate, type='prob')[,2]
	stmMask = compute_mask(sdm, climate, thresh, tol)
	stmPredictors = cbind(climate, prevalence, stmMask)		
	merge(transitionData, stmPredictors)
}


subset_transitions = function(dat, frac=0.5, max.interval=15)
{
	# subset the transition data for use in the annealing
	# first, apply the mask and drop all intervals greater than 15 years
	# then take half (or whatever fraction) of the remaining observed transitions
	# along with the same fraction of non-transitions
	intervals = dat$year2 - dat$year1
	indices = which(intervals <= max.interval & dat$sdm.mask == 1)
	transitions = with(dat[indices,], which(state1 != state2))
	notTransitions = !(transitions)
	sel = c(sample(transitions, as.integer(frac*length(transitions))),
			sample(notTransitions, as.integer(frac*length(notTransitions))))
	list(selected = dat[sel,], unselected = dat[-sel,])
}


sdm_response_curves = function(sp, dat, sdm, fname=NULL, width=10, height=7)
{
	# pull out all of the other information to just get climate vars
	otherVars = c('plot', 'year1', 'year2', 'state1', 'state2', 'lat', 'lon',
			'prevalence', 'sdm.pres', 'sdm.mask')
	climVars = colnames(dat)[!(colnames(dat) %in% otherVars)]
	
	if(is.null(fname))
	{
		dev.new(width=width, height=height)
	} else {
		pdf(file=fname, width=width, height=height)
	}
	par(mfrow=c(2,3), bty='n', oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))
	xx = matrix(0, nrow=1000, ncol=6)
	xx = as.data.frame(xx)
	colnames(xx) = climVars
	for(climVar in climVars)
	{
		xx[,climVar] = seq(-3, 3, length.out=1000)
		yy = predict(sdm, newdata=xx, type='prob')[,2]
		plot(xx[,climVar], yy, type='l', col='#016c59', xlab=climVar, 
				ylab="prob of presence", ylim=c(0,1))
		# get the SDM & STM Projection range for this var
		sdmRange = range(dat[dat$sdm.pres==1,climVar])
		projRange = range(dat[dat$sdm.mask==1,climVar])
		polygon(c(projRange, rev(projRange)), c(0,0,1,1), border=NA, col="#1c909944")
		polygon(c(sdmRange, rev(sdmRange)), c(0,0,1,1), border=NA, col="#045a8d44")
		
		xx[,climVar] = rep(0,1000)
	}
	mtext(sp, outer=T, line=0.7, cex=0.8)
	if(!is.null(fname)) dev.off()
}



plot_sdm = function(sdmDat, coords, sdm.col, legend=TRUE, add=FALSE, ...)
{
	sdmRas = make_raster(sdmDat, coords, P4S.latlon, stmMapProjection)
	plot(sdmRas, xaxt='n', yaxt='n', col=sdm.col, legend=legend, add=add, ...)
	if(!add) plot(ocean, col="white", add=TRUE)
}


plot_sdm_diagnostics = function(mapData, coord.cols, sdm.col, filename = NULL, 
		height=7, width=2*height)
{
	if(is.null(filename))
	{
		dev.new(width=width, height=height)
	} else {
		dpi = 600
		pixwidth = as.integer(dpi*width)
		pixheight = as.integer(dpi*height)
		fontsize = 15
		if(substr(filename, nchar(filename)-3, nchar(filename)) != '.png')
			filename = paste(filename, 'png', sep='.')
		png(w=pixwidth, h=pixheight, file=filename, pointsize=fontsize, res = dpi)
	}

	par(mfcol=c(1,2), mar=c(0,0,1,0), oma=c(0.5,0.5,2.5,0.5))

	# first the largest dataset, which is the mask
	plot_sdm(mapData$sdm.mask, mapData[,coord.cols], c('#ffffff00', "77777744"), legend=FALSE)
	mtext("Masks",side=3,line=0.5)

	# then the SDM presence/absence
	plot_sdm(mapData$sdm.pres, mapData[,coord.cols], c('#ffffff00', "#33333344"), 
			legend=FALSE, add=TRUE)

	# this is the actual sdm
	par(mar=c(0,0,1,3))
	plot_sdm(mapData$sdm,mapData[,coord.cols],sdm.col, zlim=c(0,1), 
			legend.args=list(side=2,text='Pr. Presence'))
	mtext("SDM",side=3,line=0.5)

	mtext(spName, side=3,line=1,outer=TRUE)

	if(!(is.null(filename)))
		dev.off()
}



#######
# MAIN LOOP ACROSS SPECIES
#######

## for(spName in speciesList)
## {
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
		
		print("  Computing ROC")
		cat("Area under the ROC curve: ", get_auc(rf.mod, sdmDat$unselected), "\n", 
			file=file.path(baseDir, 'res', 'sdm_roc.txt'))
			

		# read and merge transition dataset with the species range & SDM prob of presence
		cat("  Projecting SDM to transitions\n")
		stmData = prep_transitions(spName, baseDir, transitionClimate, rf.mod, 
				sdm.threshold, mask.tol)
		saveRDS(stmData, file.path(baseDir, 'dat', paste(spName, 'stm_data_all.rds', 
				sep='_')))
				
		# prepare the data for annealing
		annealDat = subset_transitions(stmData, frac=annealFrac, max.interval=stmMaxInterval)
		saveRDS(annealDat$selected, file.path(baseDir, 'dat', 
				paste(spName, 'stm', 'calib.rds', sep='_')))
		saveRDS(annealDat$unselected, file.path(baseDir, 'dat', 
				paste(spName, 'stm', 'valid.rds', sep='_')))


		# project SDM as well as the yes/no threshold into geographic space
		cat("  Projecting SDM geographically\n")
		sdmProjection = data.frame(lon=climGrid$lon, lat=climGrid$lat, 
				sdm=predict(rf.mod, newdata=climGrid, type='prob')[,2])
		sdmProjection = cbind(sdmProjection, compute_mask(rf.mod, climGrid, 
				sdm.threshold, mask.tol))
		saveRDS(sdmProjection, file.path(baseDir, 'res', paste(spName, 
				'sdm_grid_projection.rds', sep='_')))


		# response curve
		cat("  Producing SDM plots\n")
		sdm_response_curves(spName, stmData, rf.mod, 
				file.path(baseDir, 'img', 'sdm_response.pdf'))
				
		
		# SDM maps
		plot_sdm_diagnostics(sdmProjection, 
				which(colnames(sdmProjection) %in% c('lon', 'lat')), sdmColors, 
				file.path(baseDir, 'img', 'sdm'))
	}
}






