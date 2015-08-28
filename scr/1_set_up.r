#!/usr/bin/Rscript
# this is a quick utility script that will build the full working tree of directories
# and set up a few other things that are used globally. Should be run once on each
# computer that will be running any analyses

library(rgdal)

# create the directory tree needed by the project, if it doesn't already exist
spList = c('18032-ABI-BAL', '18037-PIN-TAE', '18086-LIR-TUL', '19049-ULM-AME', 
		'19280-QUE-NIG', '19290-QUE-ALB', '19408-QUE-RUB', '19447-QUE-VEL', 
		'19462-FAG-GRA', '19481-BET-ALL', '19489-BET-PAP', '27821-NYS-SYL', 
		'28728-ACE-RUB', '28731-ACE-SAC', '32929-FRA-PEN', '32931-FRA-AME', 
		'183295-PIC-GLA', '183302-PIC-MAR', '183385-PIN-STR', '195773-POP-TRE', 
		'19027-LIQ-STY', 'NA-CAR-ALB', '19422-QUE-STE', '19231-CAR-GLA',
		'19277-QUE-FAL', '183397-TSU-CAN', '18048-JUN-VIR', '183335-PIN-ECH',
		'505490-THU-OCC', '18034-PIC-RUB', '19051-ULM-ALA', 'NA-QUE-PRI', 
		'19242-CAR-OVA', '19288-QUE-COC', '22463-POP-GRA', '19254-JUG-NIG', 
		'19050-ULM-RUB', '23690-OXY-ARB', '19511-OST-VIR', '32945-FRA-NIG')
subDirs = c('dat', 'img', 'res')
resDirs = c('anneal', 'mcmc1', 'mcmc2', 'mcmc3')
for(sp in spList)
{
	for(subD in subDirs)
		dir.create(file.path('species', sp, subD), recursive=TRUE)
	for(subD in resDirs)
		dir.create(file.path('species', sp, 'res', subD), recursive=TRUE)
}

# directories for global data and results
for(sDir in subDirs)
	dir.create(file.path(sDir), recursive=TRUE)
	
saveRDS(spList, file.path('dat', 'speciesList.rds'))


# set up map projections
P4S.latlon = CRS("+proj=longlat +datum=WGS84")
stmMapProjection = CRS("+init=epsg:5070") # albers equal area conic NAD-83 north america
save(P4S.latlon, stmMapProjection, file="dat/map_projections.rdata")



# format and scale data
climDat = readRDS('dat/raw/plotClimate_raw.rds')
transitionClimDat = readRDS('dat/raw/transitionClimate_raw.rds')
climGrid.raw = read.csv("dat/raw/SDMClimate_grid.csv")
climGrid.raw = climGrid.raw[complete.cases(climGrid.raw),]

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
saveRDS(climScaling, "dat/climate_scaling.rds")

## add plot and year information back into the dataframes
climVars.unscaled = cbind(climDat[,c(1,2)], climDatVars)
climVars.scaled = cbind(climDat[,c(1,2)], climVars.scaled)
trClim.unscaled = cbind(transitionClimDat[,1:3], transClimDatVars)
trClim.scaled = cbind(transitionClimDat[,1:3], trClim.scaled)
climGrid.unscaled = cbind(climGrid.raw[,1:2], climGrid.unscaled)
climGrid.scaled = cbind(climGrid.raw[,1:2], climGrid.scaled)

## save climate variables
saveRDS(climVars.unscaled, "dat/plotClimate_unscaled.rds")
saveRDS(climVars.scaled, "dat/plotClimate_scaled.rds")
saveRDS(trClim.unscaled, "dat/transitionClimate_unscaled.rds")
saveRDS(trClim.scaled, "dat/transitionClimate_scaled.rds")
saveRDS(climGrid.unscaled, "dat/climateGrid_unscaled.rds")
saveRDS(climGrid.scaled, "dat/climateGrid_scaled.rds")
