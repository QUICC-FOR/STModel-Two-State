### constants
stmMaskTolerance = 1 # lat/lon tolerance for the data mask (how many degrees out of the range for projecting/calibrating
sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
		interpolate='spline', bias=1, space="rgb")(200)



prep.basedata = function()
{
	library(rgdal)
	library(raster)
	mapProj = readRDS("dat/map_projections.rds")
	ocean.up = readOGR(dsn="dat/basemaps/ne_50m_ocean", layer="ne_50m_ocean")
	ocean = spTransform(ocean.up, mapProj$projected)
	lakes = readOGR(dsn="dat/basemaps/ne_50m_lakes", layer="ne_50m_lakes")
	lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
	grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
	lakes = spTransform(lakes, mapProj$projected)
	grLakes = spTransform(grLakes, mapProj$projected)
	rivers = readOGR(dsn="dat/basemaps/ne_50m_rivers_lake_centerlines", layer="ne_50m_rivers_lake_centerlines")
	rivers = spTransform(rivers, mapProj$projected)

## 	climGrid.raw = read.csv("dat/raw/SDMClimate_grid.csv")
## 	climGrid.raw = climGrid.raw[complete.cases(climGrid.raw),]
## 	coordinates(climGrid.raw) = c('lon', 'lat')
## 	gridded(climGrid.raw) = TRUE
## 	climGrid.raw = raster(climGrid.raw)

	raster.lims = expand.grid(lon=seq(-110, -45, 0.5), lat=seq(15, 65, 0.5), val=0)
	raster.lims = rasterFromXYZ(raster.lims)
	proj4string(raster.lims) = mapProj$latlon

	set_up_baseraster = function(fname)
	{
		ras = stack(fname)
		proj4string(ras) = mapProj$latlon
		ras = crop(ras, raster.lims)
		ras = projectRaster(ras, crs=mapProj$projected)
		stack(ras)
	}
	basemap1 = set_up_baseraster('dat/basemaps/NE1_50M_SR/NE1_50M_SR.tif')
	basemap2 = set_up_baseraster('dat/basemaps/NE1_50M_SR_W/NE1_50M_SR_W.tif')
	basemap3 = set_up_baseraster('dat/basemaps/HYP_50M_SR/HYP_50M_SR.tif')
	hillshade = set_up_baseraster('dat/basemaps/GRAY_50M_SR_W/GRAY_50M_SR_W.tif')

	saveRDS(list(ocean=ocean, lakes=lakes, grLakes=grLakes, rivers=rivers,
			hillshade=hillshade, nat.earth1 = basemap1, nat.earth2 = basemap2, 
			hypso=basemap3), "dat/basedata.rds")
}
predict.stm_point = function(p, env1 = NULL, env2 = NULL)
{
	if(is.null(env1) | is.null(env2))
	{
		phi = p[1]
	} else {
		phi = p[1] + 0 * env1
		if(length(p) == 5)
			phi = phi + p[2]*env1 + p[3]*env2 + p[4]*env1^2 + p[5]*env2^2
	}
	
	result = tryCatch({
		plogis(phi)},
	error = function(e) {
		exp(phi) / (1 + exp(phi))
	})
	result
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
	lon.r = range(pres.coords[pres == 1,1]) + tol[1]*c(-1,1)
	lat.r = range(pres.coords[pres == 1,2]) + tol[2]*c(-1,1)
	inrange = function(x, r) (x >= min(r) & x <= max(r))
	mask.ind = which(inrange(new.coords[,1], lon.r) & inrange(new.coords[,2], lat.r))
	mask = rep(0, nrow(new.coords))
	mask[mask.ind] = 1
	return(mask)
}






### PLOTTING FUNCTIONS
set_up_stm_figure = function(nplots, filename=NULL, panel.width=4, panel.height=5, 
		fig.ratio=16/10, dpi=600, fontsize=15, ...)
{
	# sets up a plotting region able to hold an arbitrary mumber of plots (using mfrow)
	# great for plotting a bunch of species on one plot
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
		png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
			file=filename, pointsize=fontsize, res=dpi)
	}
	par(mfrow=c(n.panels.height,n.panels.width), ...)
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


### functions for plotting single panels
# these are intended to be called for one species to produce a single panel in a multi-
# panel plot. Thus they do no modifications of margins or anything; those should be handled
# at a higher level.

plot.sdm = function(sdmDat, sdm.col, coords=c(1,2), hill=FALSE, data.column = 3, legend=FALSE, add=FALSE, plot.ocean=TRUE, ...)
{
	# function expects a 3-column dataframe, with coordinates in first two columns
	# if this is not the case, use the coords and data.column arguments to adjust
	# this function only plots a single panel; is does no modifications of margins or 
	# anything; those should be handled at a higher level.
	require(raster)
	require(rgdal)
	if(hill | plot.ocean)
		basedata = readRDS("dat/basedata.rds")

	if(hill)
	{
		plotRGB(basedata$hillshade)
		add = TRUE		
	}
	if(!(identical(coords, c(1,2))))
		sdmDat = sdmDat[, c(coords, data.column)]
	sdmRas = rasterFromXYZ(sdmDat)
	plot(sdmRas, xaxt='n', yaxt='n', col=sdm.col, legend=legend, add=add, ...)

	if(plot.ocean)
	{
		plot(basedata$ocean, col="white", add=TRUE, lwd=0.5)
		plot(basedata$grLakes, col="white", add=TRUE, lwd=0.5)
	}
}

