### constants
stmMaskTolerance = 1 # lat/lon tolerance for the data mask (how many degrees out of the range for projecting/calibrating
sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
		interpolate='spline', bias=1, space="rgb")(200)

### data
library(rgdal)
load("dat/map_projections.rdata")
ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)


### general functions
# deprecated - use the predict.stm_point function instead
## compute_e = function(p, env1, env2) plogis(p[6] + env1*p[7] + env2*p[8] + env1^2*p[9] + 
## 		env2^2*p[10])
## compute_c = function(p, env1, env2) plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + 
## 		env2^2*p[5])
predict.stm_point = function(p, env1 = NA, env2 = NA)
{
	phi = p[1]
	if(length(p == 5)) phi = phi + p[2]*env1 + p[3]*env2 + p[4]*env1^2 + p[5]*env2^2
	plogis(phi)
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

plot_sdm = function(sdmDat, coords, sdm.col, legend=FALSE, add=FALSE, plot.ocean=TRUE, ...)
{
	sdmRas = make_raster(sdmDat, coords, P4S.latlon, stmMapProjection)
	plot(sdmRas, xaxt='n', yaxt='n', col=sdm.col, legend=legend, add=add, ...)
	if(plot.ocean) plot(ocean, col="white", add=TRUE)
}


