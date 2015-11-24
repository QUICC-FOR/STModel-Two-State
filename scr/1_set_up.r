#!/usr/bin/env Rscript
# set up data files etc

## depends
##   dat/raw

## produces
##   dat/clim/climate_scaling.rds

dir.create(file.path('dat', 'clim'), recursive=TRUE)


## species names
# create the directory tree needed by the project, if it doesn't already exist
spList = c(
# temperate species already done
'28728-ACE-RUB', '28731-ACE-SAC', '19462-FAG-GRA', '32931-FRA-AME',
'32945-FRA-NIG', '19287-QUE-MAC', '19408-QUE-RUB', '183397-TSU-CAN',

# transitional species already done
'19481-BET-ALL', '183375-PIN-RES', '183385-PIN-STR', '22463-POP-GRA',

# boreal species already done
'18032-ABI-BAL', '19489-BET-PAP', '183412-LAR-LAR', '183295-PIC-GLA', '183302-PIC-MAR', 
'18034-PIC-RUB', '183319-PIN-BAN', '195773-POP-TRE', '505490-THU-OCC'
)
saveRDS(spList, file.path('dat', 'speciesList.rds'))


# set up map projections
if(!file.exists("dat/map_projections.rdata"))
{
	tryCatch(
	{
		library(rgdal)
		P4S.latlon = CRS("+proj=longlat +datum=WGS84")
		stmMapProjection = CRS("+init=epsg:5070") # albers equal area conic NAD-83 north america
		saveRDS(list(latlon=P4S.latlon, projected=stmMapProjection, 
				proj.name="albers equal area conic NAD-83 north america"), 
				file="dat/map_projections.rds")
	}, error=function(e)
	{
		msg = paste("Could not save map projections to dat/map_projections.rdata\n",
			"error:\n", e)
		warning(msg)
	})
}



# format and scale data
climDat = readRDS('dat/raw/plotClimate_raw.rds')
transitionClimDat = readRDS('dat/raw/transitionClimate_raw.rds')
climGrid.raw = read.csv("dat/raw/SDMClimate_grid.csv")


# remove some empty spaces in Ontario and Newfoundland from the rasters
# this was done because there is no plot coverage there and there were inappropriate
# projections happening
remove_black_holes = function(x)
{	
	require(raster)
	require(rgdal)
	require(sp)
	
	mapProj = readRDS("dat/map_projections.rds")

	x = x[complete.cases(x),]
	coordinates(x) = x[,c('lon', 'lat')]
	gridded(x) = TRUE
	x = brick(x)
	proj4string(x) = mapProj$latlon
	
	# remove great lakes from the projections	
	lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
	lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
	grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
	grLakes.ras = rasterize(grLakes,x)
	x[!is.na(grLakes.ras)] = NA

	# ontario
	e.lim = -79.4
	w.lim = -100
	n.lim = 54
	s.lim = 50.2
	# newfoundland
	nf.e.lim = -53
	nf.w.lim = -59.5
	nf.n.lim = 50.3
	nf.s.lim = 47
	#florida keys
	fk.w = -79.5
	fk.n = 30
	# random atlantic islands
	i1.w=-68
	i1.n=35
	i2.w=-61
	i2.n=44.5
	
	x[x$lon < e.lim & x$lon > w.lim & x$lat > s.lim & x$lat < n.lim] = NA
	x[x$lon > nf.w.lim & x$lon < nf.e.lim & x$lat > nf.s.lim & x$lat < nf.n.lim] = NA
	x[x$lon > fk.w & x$lat < fk.n] = NA
	x[x$lon > i1.w & x$lat < i1.n] = NA
	x[x$lon > i2.w & x$lat < i2.n] = NA
	x = as.data.frame(x)
	x[complete.cases(x),]
	
}
climGrid.raw = remove_black_holes(climGrid.raw)

# drop the variable that is causing problems (gdd period 2)
# inspection of the data suggested that this variable is corrupted; it is NOT GDD
climDat = climDat[,-7]

# PCA
## var.pca = dudi.pca(climDat, scannf=FALSE, nf = 5)
## print("PCA cumulative variance explained:")
## print(var.pca$eig / sum(var.pca$eig))
## varCor = cor(climDat[-c(1,2)])
## contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

# procedure:
# select variables that are relatively uncorrelated to each other
# use PCA to find variables that explain unique variance
# start with mean annual pp and temp (as overall representative)
# in practice, there are 3 uncorrelated temp and 4 precip variables
# we drop one precip variable to have 3 of each

# it is commented out now, because it is really an interactive procedure; provided here
# for documentation
## nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),
## 		"annual_mean_temp"]), names(varCor[which(abs(varCor[,"tot_annual_pp"])<0.7),
## 		"tot_annual_pp"]))
## print(contrib[nonCorVars,])

selectedVars = c('annual_mean_temp', 'mean_diurnal_range', 'mean_temp_wettest_quarter',
		'tot_annual_pp', 'pp_seasonality', 'pp_warmest_quarter')
		
# extract only the climate variables we need
climCols = which(colnames(climDat) %in% selectedVars)
climDatVars = climDat[, climCols]
trClimCols = which(colnames(transitionClimDat) %in% selectedVars)
transClimDatVars = transitionClimDat[, trClimCols]
cgDatCols = which(colnames(climGrid.raw) %in% selectedVars)
climGrid.unscaled = climGrid.raw[,cgDatCols]

## scale the variables and save the scaling
climVars.scaled = scale(climDatVars)
trClim.scaled = scale(transClimDatVars, center = attr(climVars.scaled, "scaled:center"),
		scale = attr(climVars.scaled, "scaled:scale"))
climGrid.scaled = scale(climGrid.unscaled, center = attr(climVars.scaled, "scaled:center"),
		scale = attr(climVars.scaled, "scaled:scale"))
climScaling = list(center = attr(climVars.scaled, "scaled:center"),
		scale = attr(climVars.scaled, "scaled:scale"))
saveRDS(climScaling, "dat/clim/climate_scaling.rds")

## add plot and year information back into the dataframes
climVars.unscaled = cbind(climDat[,c(1,2)], climDatVars)
climVars.scaled = cbind(climDat[,c(1,2)], climVars.scaled)
trClim.unscaled = cbind(transitionClimDat[,1:3], transClimDatVars)
trClim.scaled = cbind(transitionClimDat[,1:3], trClim.scaled)
climGrid.unscaled = cbind(climGrid.raw[,1:2], climGrid.unscaled)
climGrid.scaled = cbind(climGrid.raw[,1:2], climGrid.scaled)

## project everything
## for map data, need to project it as a raster then convert back to df; this ensures
## that it stays as a gridded data type for mapping later 
## (faster to manipulate as df than as raster)
df.project = function(r, coords = c('lon', 'lat'))
{
	require(raster)
	require(rgdal)
	require(sp)
	coordinates(r) = coords
	gridded(r) = TRUE
	r = brick(r)
	proj4string(r) = P4S.latlon
	r = projectRaster(r, crs=stmMapProjection)
	rdf = as.data.frame(r)
	rdf = cbind(rdf, coordinates(r))
	rdf[complete.cases(rdf),]
}
climGrid.scaled.pr = df.project(climGrid.scaled)
climGrid.unscaled.pr = df.project(climGrid.unscaled)
plot.locs = readRDS('dat/raw/plotLocations.rds')
if(!('x' %in% colnames(plot.locs)) | !('y' %in% colnames(plot.locs)))
{
	require(sp)
	coordinates(plot.locs) = plot.locs[,c('lon', 'lat')]
	proj4string(plot.locs) = P4S.latlon
	plot.locs = spTransform(plot.locs, stmMapProjection)
	plot.locs = as.data.frame(plot.locs)
	colnames(plot.locs)[4:5] = c('x', 'y')
	saveRDS(plot.locs, "dat/raw/plotLocations.rds")
}

## save climate variables
saveRDS(climVars.unscaled, "dat/clim/plotClimate_unscaled.rds")
saveRDS(climVars.scaled, "dat/clim/plotClimate_scaled.rds")
saveRDS(trClim.unscaled, "dat/clim/transitionClimate_unscaled.rds")
saveRDS(trClim.scaled, "dat/clim/transitionClimate_scaled.rds")
saveRDS(climGrid.unscaled.pr, "dat/clim/climateGrid_unscaled.rds")
saveRDS(climGrid.unscaled, "dat/clim/climateGrid_unscaled_unprojected.rds")
saveRDS(climGrid.scaled.pr, "dat/clim/climateGrid_scaled.rds")
