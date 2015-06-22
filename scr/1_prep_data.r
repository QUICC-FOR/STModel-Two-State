#!/usr/bin/Rscript
library(ade4)
library(raster)


# set up map projections
P4S.latlon = CRS("+proj=longlat +datum=WGS84")
stmMapProjection = CRS("+init=epsg:6350") # albers equal area conic NAD-83 north america
save(P4S.latlon, stmMapProjection, file="dat/map_projections.rdata")


climDat = readRDS('dat/raw/plotClimate_raw.rds')
transitionClimDat = readRDS('dat/raw/transitionClimate_raw.rds')
climGrid.raw = read.csv("dat/raw/SDMClimate_grid.csv")
climGrid.raw = climGrid.raw[complete.cases(climGrid.raw),]

# drop the variable that is causing problems (gdd period 2)
# inspection of the data suggested that this variable is corrupted; it is NOT GDD
climDat = climDat[,-7]

# PCA
var.pca = dudi.pca(climDat, scannf=FALSE, nf = 5)
print("PCA cumulative variance explained:")
print(var.pca$eig / sum(var.pca$eig))
varCor = cor(climDat[-c(1,2)])
contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

# procedure:
# select variables that are relatively uncorrelated to each other
# use PCA to find variables that explain unique variance
# start with mean annual pp and temp (as overall representative)
# in practice, there are 3 uncorrelated temp and 4 precip variables
# we drop one precip variable to have 3 of each
nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),
		"annual_mean_temp"]), names(varCor[which(abs(varCor[,"tot_annual_pp"])<0.7),
		"tot_annual_pp"]))
print(contrib[nonCorVars,])

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
saveRDS(climScaling, "dat/climate_scaling.rds")

## add plot and year information back into the dataframes
climVars.scaled = cbind(climDat[,c(1,2)], climVars.scaled)
trClim.scaled = cbind(transitionClimDat[,1:3], trClim.scaled)
climGrid.unscaled = cbind(climGrid.raw[,1:2], climGrid.unscaled)
climGrid.scaled = cbind(climGrid.raw[,1:2], climGrid.scaled)

# make the climGrids into projected raster brick objects
make_brick = function(x)
{
	coordinates(x) = c('lon', 'lat')
	gridded(x) = TRUE
	x.ras = brick(x)
	setMinMax(x.ras)
	proj4string(x.ras) = P4S.latlon
	projectRaster(x.ras, crs=stmMapProjection)
}
climGrid.unscaled = make_brick(climGrid.unscaled)
climGrid.scaled = make_brick(climGrid.scaled)


## save climate variables
saveRDS(climVars.scaled, "dat/plotClimate_scaled.rds")
saveRDS(trClim.scaled, "dat/transitionClimate_scaled.rds")
saveRDS(climGrid.unscaled, "dat/climateGrid_unscaled.rds")
saveRDS(climGrid.scaled, "dat/climateGrid_scaled.rds")
