### general functions

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






### PLOTTING FUNCTIONS

### functions for plotting single panels
# these are intended to be called for one species to produce a single panel in a multi-
# panel plot. Thus they do no modifications of margins or anything; those should be handled
# at a higher level.

