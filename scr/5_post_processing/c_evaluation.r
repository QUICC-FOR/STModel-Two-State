#!/usr/bin/env Rscript

# prepares temp plot data and evaluates the SDM/STM for all species specified
# gets posterior tss/roc at the map scale

# depends:
#    dat/clim/climate_scaling.rds
#    dat/raw/temp_plots/tp_climData_reshaped.csv
#    dat/raw/temp_plots/tp_plotinfoData.csv
#    dat/raw/temp_plots/tp_treeData_allSpecies.csv
#    dat/pp_grid_tall.rds
#    dat/tp_grid_tall.rds
#    dat/transition/*
#    dat/presence/*
#    res/sdm/*
#    res/posterior/*

# produces:
#    dat/tempPlot_presence.rds
#    res/eval/sdm_eval.rds
#    res/eval/stm_eval.rds

library(doParallel)
library(foreach)
library(randomForest)
library(biomod2)
numCores = detectCores()
registerDoParallel(cores=numCores)
modName = '0'
posterior.n = 1000

suppressWarnings(
{
	dir.create(file.path('res', 'eval'), recursive=TRUE)
	dir.create(file.path('dat', 'eval'), recursive=TRUE)
})

speciesList = readRDS("dat/speciesList.rds")

source('scr/stm_functions.r')

arg = commandArgs(TRUE)
if(length(arg) > 0)
{
	sps = which(speciesList %in% arg)
	if(length(sps) == 0) 
	{
		warning("No species specified on command line; falling back to default")
	} else 
		speciesList = speciesList[which(speciesList %in% arg)]
}

get_temp_plots = function()
{
	require(reshape2)
	if(file.exists("dat/tempPlot_presence.rds"))
	{
		tempPlotsAll = readRDS("dat/tempPlot_presence.rds")
	} else {
		# process temp plot climate
		climScale = readRDS("dat/clim/climate_scaling.rds")
		tp_clim = read.csv("dat/raw/temp_plots/tp_climData_reshaped.csv", dec=',', sep=';')
		colsToKeep = c('plot_id', 'year_measured')
		tp_clim = cbind(tp_clim[,colsToKeep], scale(tp_clim[,names(climScale$center)], 
				center=climScale$center, scale=climScale$scale))
		tp_clim = tp_clim[complete.cases(tp_clim),]
		tempPlots = read.csv("dat/raw/temp_plots/tp_plotinfoData.csv", dec='.', sep=';')
		tp_species = read.csv("dat/raw/temp_plots/tp_treeData_allSpecies.csv", dec='.', sep=';')
		treeDat = merge(tp_species, tempPlots, all.x=TRUE)

		# drop unneeded stuff and save
		treeDat = treeDat[-(which(treeDat$id_spe == "")),-7]
		sampleDat = dcast(treeDat, plot_id + year_measured + latitude + longitude ~ id_spe, fill = 0, 
				value.var = "basal_area", fun.aggregate = function(x) as.integer(sum(x) > 0))
		tempPlotsAll = merge(sampleDat, tp_clim, all = 'T', by=c("plot_id", "year_measured"))
		tempPlotsAll = tempPlotsAll[complete.cases(tempPlotsAll),]
		cn = colnames(tempPlotsAll)
		
		# project the points
		mapProj = readRDS("/Users/mtalluto/Dropbox/work/projects_git/STModel-Two-State/dat/map_projections.rds")
		coordinates(tempPlotsAll) = tempPlotsAll[,c('longitude', 'latitude')]
		proj4string(tempPlotsAll) = mapProj$latlon
		tempPlotsAll = spTransform(tempPlotsAll, mapProj$projected)
		tempPlotsAll = as.data.frame(tempPlotsAll)
		names(tempPlotsAll) = c(cn, 'x', 'y')
		saveRDS(tempPlotsAll, "dat/tempPlot_presence.rds")
	}
	tempPlotsAll
}



######
# process annual climate rasters
make.brick.list = function(dat)
{
	require(raster)
	require(rgdal)
	require(reshape2)
	
	ras.wide = dcast(dat, year_measured + lon + lat ~ biovar, value.var = 'val')
	ras.wide = ras.wide[complete.cases(ras.wide),]
	
	# for some reason these were multiplied by 10
	ras.wide[,'mean_diurnal_range'] = ras.wide[,'mean_diurnal_range'] / 10
	ras.wide[,'mean_temp_wettest_quarter'] = ras.wide[,'mean_temp_wettest_quarter'] / 10

	# scale the climate data
	climScale = readRDS("dat/clim/climate_scaling.rds")
	ras.wide.scale = ras.wide
	for(v in names(climScale$center))
		ras.wide.scale[,v] = scale(ras.wide.scale[,v], center = climScale$center[v], scale=climScale$scale[v])

	# make a list by year, each element is a projected raster brick for a given year; 
	# each layer is a biovar
	foreach(year=unique(ras.wide.scale$year_measured), .inorder=TRUE, 
			.final=function(x) { 
				names(x) = unique(ras.wide.scale$year_measured)
				x
			}) %dopar%
	{
		dfsub = ras.wide.scale[ras.wide.scale$year_measured == year,2:ncol(ras.wide.scale)]
		coordinates(dfsub) = c('lon', 'lat')
		gridded(dfsub) = TRUE
		climBrick = brick(dfsub)
		mapProj = readRDS("dat/map_projections.rds")
		proj4string(climBrick) = mapProj$latlon
		projectRaster(climBrick, crs=mapProj$projected)
	}
}

permPlotRasterFile = 'dat/clim/permPlot_annual_clim_raster.rds'
tempPlotRasterFile = 'dat/clim/tempPlot_annual_clim_raster.rds'
if(!file.exists(permPlotRasterFile))
{
	permPlot.annualClimRaster = make.brick.list(readRDS('dat/raw/permPlot_grid_tall.rds'))
	saveRDS(permPlot.annualClimRaster, permPlotRasterFile)
} else {
	permPlot.annualClimRaster = readRDS(permPlotRasterFile)
}
if(!file.exists(tempPlotRasterFile))
{
	tempPlot.annualClimRaster = make.brick.list(readRDS('dat/raw/tempPlot_grid_tall.rds'))
	saveRDS(tempPlot.annualClimRaster, tempPlotRasterFile)
} else {
	tempPlot.annualClimRaster = readRDS(tempPlotRasterFile)
}


########
#
#   Process individual species
#
#######
plotLocs = readRDS('dat/raw/plotLocations.rds')
rfvars = c('annual_mean_temp', 'mean_diurnal_range', 'mean_temp_wettest_quarter', 
		'pp_seasonality', 'pp_warmest_quarter', 'tot_annual_pp')

intersect_raster_cells = function(plots, rasBrick, spName)
{
	require(raster)
	do.call(rbind,
		lapply(unique(plots$year), function(year) {
			pl = plots[plots$year == year,c(spName, 'x', 'y')]		
			coordinates(pl) = c('x', 'y')
			rb = rasBrick[[as.character(year)]]
			rb$obs = rasterize(pl, rb, 
					fun=function(x, ...) as.integer(sum(x, ...) > 0))[[2]]
			res = cbind(getValues(rb), coordinates(rb))
			res[complete.cases(res),]
		})
	)
}


tempPlotsAll = get_temp_plots()

# evaluate the SDM and create a list of raster cells for validation (at the grid/map scale)
eval.rf = mclapply(speciesList, function(spName)
## for(spName in speciesList)
{
	# read in calibration data and find points that were only sampled once
	# (and thus not used to fit the STM)
	trans = readRDS(file.path('dat', 'transition', paste0(spName, '_transitions.rds')))
	pres = readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))
	singlePlots = pres[which(!(pres$plot_id %in% trans$plot)),]
	singlePlots = merge(singlePlots, plotLocs)
	singlePlots = singlePlots[complete.cases(singlePlots),]

	# set up temp plots
	tempPlots = tempPlotsAll[,c('plot_id', 'year_measured', spName, 'x', 'y', 'longitude', 'latitude')]
	tempPlots = tempPlots[complete.cases(tempPlots),]

	# for each year, compute obs and converte to data frame with complete cases only
 	mapValidCells = rbind(
		intersect_raster_cells(tempPlots, tempPlot.annualClimRaster, spName),
		intersect_raster_cells(singlePlots, permPlot.annualClimRaster, spName)
	)
	mapValidCells = as.data.frame(mapValidCells)
	# get SDM fit
	RFMod = readRDS(file.path('res', 'sdm', paste0(spName, '_rf_sdm.rds')))
	mapValidCells$fit.rf = predict(RFMod, newdata = mapValidCells[,rfvars], type='prob')[,2]
	saveRDS(mapValidCells, file.path('dat', 'eval', paste0(spName, '_evalCells.rds')))
	
	# compute tss/auc and return
	list(validation.data = mapValidCells, validation.results =  with(mapValidCells, 
			rbind(Find.Optim.Stat(Stat="TSS", Fit=(1000*fit.rf), Obs=obs),
			Find.Optim.Stat(Stat="ROC", Fit=(1000*fit.rf), Obs=obs)))
	)
}, mc.cores=numCores)
## }
names(eval.rf) = speciesList
cat("sdm evaluation done\n")


#### now evaluate the STM

stm.eval = foreach(spName=speciesList, .final = function(x) {names(x) = speciesList; x}) %do%
{
	
	# get STM fit
	samples = readRDS(file.path('res', 'posterior', paste(spName, modName, 'samples.rds', sep='_')))
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]
	mapValidCells = eval.rf[[spName]][['validation.data']]
	env1 = mapValidCells$annual_mean_temp
	env2 = mapValidCells$tot_annual_pp
	
	stm.reps = foreach(pars = iter(samples, by='row'), .packages='biomod2', 
			.combine=rbind) %dopar%
	{
		if(length(pars) == 2)
		{
			cp = pars[1]
			ep = pars[2]
		} else {
			cp = pars[1:5]
			ep = pars[6:10]
		}
		C = predict.stm_point(cp, env1, env2)
		E = predict.stm_point(ep, env1, env2)
		as.integer(C > E)
	}
	stm.fit = colSums(stm.reps) / ncol(stm.reps)
	tss = Find.Optim.Stat(Stat="TSS", Fit=(1000*stm.fit), Obs=mapValidCells$obs)
	roc = Find.Optim.Stat(Stat="ROC", Fit=(1000*stm.fit), Obs=mapValidCells$obs)
	c(tss=tss[1], tss.thresh=tss[2]/1000, roc=roc[1], roc.thresh=roc[2]/1000)
}
for(spName in speciesList) {
	saveRDS(stm.eval[[spName]], file.path('res', 'eval', paste0(spName, '_stm_eval.rds')))
	saveRDS(eval.rf[[spName]], file.path('res', 'eval', paste0(spName, '_sdm_eval.rds')))
}
cat("stm evaluation done\n")

