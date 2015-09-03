#!/usr/bin/Rscript
library(randomForest)
library(raster)
library(rgdal)

# some constants
sdm.threshold = 0.25
annealFrac = 0.5
stmMaxInterval = 15
mask.tol = 2 # lat/lon tolerance for the data mask (how many degrees out of the range for projecting/calibrating
sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
		interpolate='spline', bias=1, space="rgb")(200)


speciesList = readRDS('dat/speciesList.rds')
transitionClimate = readRDS('dat/transitionClimate_scaled.rds')
climGrid = readRDS('dat/climateGrid_scaled.rds')
plotLocs = readRDS('dat/plotLocations.rds')
transitionClimate = merge(transitionClimate, plotLocs, by.x=1, by.y=1)
load('dat/map_projections.rdata')
ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)



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



sdm_pres = function(sdm, newdata, threshold)
{
	# return a presence absence vector using the predictions of an sdm and a threshold
	predictions = predict(sdm, newdata = newdata, type='prob')[,2]
	as.integer(predictions >= threshold)

}

stm_mask = function(new.coords, pres, pres.coords, tol = c(10,10))
{
	# computes a data selection mask based on the lat/long of observed presences
	# new.coords: the coordinates to apply the mask to
	# pres: a presence absence vector
	# pres.coords: data frame or matrix with x,y coords of the pres data
	# tol: tolerance in x,y degrees
	# returns a mask vector of same length as newdata showing the data limits within 
	# tol degrees of lon/lat
	# if tol is length 1, it will be applied in both dimensions
	# if length 2, then the first will be longitude, second latitude

	if(length(tol) == 1) tol = rep(tol,2)
	lon.r = range(pres.coords[sdm.pres == 1,1]) + tol[1]*c(-1,1)
	lat.r = range(pres.coords[sdm.pres == 1,2]) + tol[2]*c(-1,1)
	inrange = function(x, r) (x >= min(r) & x <= max(r))
	mask.ind = which(inrange(new.coords[,1], lon.r) & inrange(new.coords[,2], lat.r))
	mask = rep(0, nrow(new.coords))
	mask[mask.ind] = 1
	return(mask)
}


prep_transitions = function(sp, baseDir, climate, presence, sdm, thresh, tol, 
		plot.locations, coord.names=c('lon','lat'))
{
	by.tr = 'plot'
	by.loc = 'plot_id'
	transitionData = readRDS(file.path(baseDir, 'dat', paste(sp, 'transitions.rds', sep='_')))
	transitionData = merge(transitionData, plot.locations, by.x=by.tr, by.y=by.loc)
	prevalence = predict(sdm, newdata = climate, type='prob')[,2]
	sdmPres = sdm_pres(sdm, climate, thresh)
	stmMask = stm_mask(climate[,coord.names], presence[,sp], presence[,coord.names], tol)
	stmPredictors = cbind(climate, prevalence, sdmPres, stmMask)
	colnames(stmPredictors) = c(colnames(climate), 'prevalence', 'sdm.presence', 'stm.mask')
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



plot_sdm = function(sdmDat, coords, sdm.col, legend=FALSE, add=FALSE, plot.ocean=TRUE, ...)
{
	sdmRas = make_raster(sdmDat, coords, P4S.latlon, stmMapProjection)
	plot(sdmRas, xaxt='n', yaxt='n', col=sdm.col, legend=legend, add=add, ...)
	if(plot.ocean) plot(ocean, col="white", add=TRUE)
}



#######
# MAIN LOOPS ACROSS SPECIES
#   3 passes:
#      1: read and prep data, make response curves
#      2: make a big SDM figure for all spp
#      3: make a big mask figure for all spp
#######

presences = list()
transitions = list()
projections = list()

for(spName in speciesList)
{
	cat(paste("Starting species", spName, '\n'))

	baseDir = file.path('species', spName)
	rfModFilename = file.path(baseDir, 'res', paste(spName, 'sdm.rds', sep='_'))
	rf.mod = readRDS(rfModFilename)
	
	# read and merge transition dataset with the species range & SDM prob of presence
	cat("  Projecting SDM to transitions\n")
	presDat = readRDS(file.path(baseDir, 'dat', paste(spName, 'presence.rds', sep='_')))
	presDat = merge(presDat, plotLocs)
	prLocs = presDat[presDat[,spName] == 1, c('lon', 'lat', spName)]
	coordinates(prLocs) = c('lon', 'lat')
	presences[[spName]] = prLocs
	stmData = prep_transitions(spName, baseDir, transitionClimate, presDat, rf.mod, 
				sdm.threshold, mask.tol, plotLocs)
	saveRDS(stmData, file.path(baseDir, 'dat', paste(spName, 'stm_data_all.rds', 
			sep='_')))

	# prepare the data for annealing
	annealDat = subset_transitions(stmData, frac=annealFrac, max.interval=stmMaxInterval)
	saveRDS(annealDat$selected, file.path(baseDir, 'dat', 
			paste(spName, 'stm', 'calib.rds', sep='_')))
	saveRDS(annealDat$unselected, file.path(baseDir, 'dat', 
			paste(spName, 'stm', 'valid.rds', sep='_')))

	# massage data for mapping
	class.trans = function(st1, st2)
	{
		res = character(length(st1))
		res[st1 == 0 & st2 == 0] = 'absence'
		res[st1 == 0 & st2 == 1] = 'colonization'
		res[st1 == 1 & st2 == 1] = 'presence'
		res[st1 == 1 & st2 == 0] = 'extinction'
		res
	}
	transLocs = data.frame(
		transition = factor(c(with(annealDat$selected, class.trans(state1, state2)), 
			with(annealDat$unselected, class.trans(state1, state2)))),
		lon = c(annealDat$selected$lon, annealDat$unselected$lon),
		lat = c(annealDat$selected$lat, annealDat$unselected$lat))
	coordinates(transLocs) = c('lon', 'lat')
	transitions[[spName]] = transLocs
	
	# project SDM as well as the yes/no threshold into geographic space
	cat("  Projecting SDM geographically\n")
	sdmProjection = data.frame(lon=climGrid$lon, lat=climGrid$lat, 
			sdm=predict(rf.mod, newdata=climGrid, type='prob')[,2])
	sdmProjection$sdm.pres = sdm_pres(rf.mod, climGrid, sdm.threshold)
	sdmProjetion$stm.mask = with(sdmProjection, stm_mask(cbind(lon, lat), presDat[,sp], 
			presDat[,c('lon', 'lat')], mask.tol))
	saveRDS(sdmProjection, file.path(baseDir, 'res', paste(spName, 
			'sdm_grid_projection.rds', sep='_')))
	projections[[spName]] = sdmProjection


	# response curves
	sdm_response_curves(spName, stmData, rf.mod, file.path(baseDir, 'img', 'sdm_response.pdf'))
}


# now, set up plotting region (png, par, etc)
# then loop across species
# then plot the sdm for each
# try plotting presences on the map for kicks

# then set up another plotting region
# then loop across species
# then plot the mask for each species
# also plot the transitions on the map
# last panel will be empty, but with a legend


set_up_figure = function(nplots, filename=NULL, panel.width=4, panel.height=5, 
		fig.ratio=16/10, dpi=600, fontsize=15, ...)
{
	n.panels.height = round(sqrt(nplots/fig.ratio),0)
	n.panels.width = ceiling(nplots/n.panels.height)
	figure.width = panel.width*n.panels.width
	figure.height = panel.height*n.panels.height

	if(is.null(filename) || nchar(filename) == 0)
	{
		dev.new(width=figure.width, height=figure.height)
	} else {
		fontsize = 15
		if(substr(filename, nchar(filename)-3, nchar(filename)) != '.png')
			filename = paste(filename, 'png', sep='.')
		png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height)
			file=filename, pointsize=fontsize, res=dpi)
	}
	par(mfrow=c(n.panels.height,n.panels.height), ...)
}

clean_up_figure = function(filename)
{
	if(is.null(filename) || nchar(filename) == 0)
	{
		return()
	} else {
		dev.off()
	}
}


# set up the SDM plot
coord.cols=c('lon', 'lat')
fname = file.path('img', 'sdms.png')
set_up_figure(length(speciesList) + 1, fname, mar=c(0,0,1,0), oma=c(0.5,0.5,2.5,0.5))
for(spName in speciesList)
{
	mapData = projections[[spName]]
	plot_sdm(mapData$sdm, mapData[,coord.cols], sdmColors, zlim=c(0,1))
	mtext(spName,side=3,line=0.5)
}
# plot a legend
plot_sdm(mapData$sdm, mapData[,coord.cols], sdmColors, zlim=c(0,1), legend=TRUE, 
		legend.only=TRUE, plot.ocean=FALSE)
clean_up_figure(fname)


# now plot the masks
fname = file.path('img', 'masks.png')
trans.colors = c('#1b9e77','#d95f02','#7570b3','#e7298a')
set_up_figure(length(speciesList) + 1, fname, mar=c(0,0,1,0), oma=c(0.5,0.5,2.5,0.5))
for(spName in speciesList)
{
	mapData = projections[[spName]]
	# first the largest dataset, which is the mask
	plot_sdm(mapData$sdm.mask, mapData[,coord.cols], c('#ffffff00', "77777744"), 
			legend=FALSE, plot.ocean=FALSE)
	# then the SDM presence/absence
	plot_sdm(mapData$sdm.pres, mapData[,coord.cols], c('#ffffff00', "#33333344"), 
			legend=FALSE, add=TRUE)
	mtext("Masks",side=3,line=0.5)
	points(transitions, col='black', bg=trans.cols, pch=21, cex=0.4)
}
plot(0,0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1) ylim=c(0,1))
legend(0,1, legend=c('sdm presence', 'data mask', levels(transitions[[spName]])), 
		pt.bg=c('#777777', '#333333', trans.colors), col='black', cex=2, 
		pch=c(22, 22, rep(21, length(levels(transitions[[spName]])))))
clean_up_figure(fname)


