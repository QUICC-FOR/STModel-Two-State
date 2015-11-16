#!/usr/bin/env Rscript

# prepares temp plot data and evaluates the SDM for all species
# depends:
#    dat/climate_scaling.rds
#    dat/raw/temp_plots/tp_climData_reshaped.csv
#    dat/raw/temp_plots/tp_plotinfoData.csv
#    dat/raw/temp_plots/tp_treeData_allSpecies.csv
#    dat/pp_grid_tall.rds
#    dat/transition/*
#    dat/presence/*
#    res/sdm/*

# produces:
#    dat/tempPlot_presence.rds
#    dat/eval/*
#    res/eval/sdm_eval.rds

library(raster)
library(reshape2)
library(foreach)
library(doParallel)
library(randomForest)
library(parallel)
library(biomod2)

numCores = 7
registerDoParallel(cores=numCores)

suppressWarnings(
{
	dir.create(file.path('dat', 'eval'), recursive=TRUE)
	dir.create(file.path('res', 'eval'), recursive=TRUE)
})

climScale = readRDS("dat/climate_scaling.rds")
speciesList = readRDS("dat/speciesList.rds")

# process temp plot climate
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
saveRDS(tempPlotsAll, "dat/tempPlot_presence.rds")



######
# process annual climate grids
make_grid_list = function(gr)
{
	ras.wide = dcast(gr, year_measured + lon + lat ~ biovar, value.var = 'val')
	ras.wide = ras.wide[complete.cases(ras.wide),]
	ras.wide[,'mean_diurnal_range'] = ras.wide[,'mean_diurnal_range'] / 10
	ras.wide[,'mean_temp_wettest_quarter'] = ras.wide[,'mean_temp_wettest_quarter'] / 10
	ras.wide.scale = ras.wide
	for(v in names(climScale$center))
		ras.wide.scale[,v] = scale(ras.wide.scale[,v], center = climScale$center[v], scale=climScale$scale[v])

	# make a list by year, each element is a raster brick for a given year; each layer is a biovar
	foreach(year=unique(ras.wide.scale$year_measured), .inorder=TRUE, 
			.final=function(x) { 
				names(x) = unique(ras.wide.scale$year_measured)
				x
			}) %dopar%
	{
		dfsub = ras.wide.scale[ras.wide.scale$year_measured == year,2:ncol(ras.wide.scale)]
		coordinates(dfsub) = 1:2
		gridded(dfsub) = TRUE
		brick(dfsub)
	}
}
ppg = readRDS('dat/pp_grid_tall.rds')
tpg = readRDS('dat/tp_grid_tall.rds')

ppgClimList = make_grid_list(ppg)
tpgClimList = make_grid_list(tpg)

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
	do.call(rbind,
		lapply(unique(plots$year), function(year) {
			pl = plots[plots$year == year,c(spName, 'lon', 'lat')]		
			coordinates(pl) = c('lon', 'lat')
			rb = rasBrick[[as.character(year)]]
			rb$obs = rasterize(pl, rb, 
					fun=function(x, ...) as.integer(sum(x, ...) > 0))[[2]]
			res = getValues(rb)
			res[complete.cases(res),]
		})
	)
}

eval.rf = mclapply(speciesList, function(spName)
{
	# read in calibration data and find points that were only sampled once
	# (and thus not used to fit the STM)
	trans = readRDS(file.path('dat', 'transition', paste0(spName, '_transitions.rds')))
	pres = readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))
	singlePlots = pres[which(!(pres$plot_id %in% trans$plot)),]
	singlePlots = merge(singlePlots, plotLocs)
	singlePlots = singlePlots[complete.cases(singlePlots),]

	# set up temp plots
	tempPlots = tempPlotsAll[,c('plot_id', 'year_measured', spName, 'latitude', 'longitude')]
	names(tempPlots) = names(singlePlots)
	tempPlots = tempPlots[complete.cases(tempPlots),]

	# for each year, compute obs and converte to data frame with complete cases only
 	mapValidCells = rbind(
		intersect_raster_cells(tempPlots, tpgClimList, spName),
		intersect_raster_cells(singlePlots, ppgClimList, spName)
	)
	mapValidCells = as.data.frame(mapValidCells)
	# get SDM fit
	RFMod = readRDS(file.path('res', 'sdm', paste0(spName, '_rf_sdm.rds')))
	mapValidCells$fit.rf = predict(RFMod, newdata = mapValidCells[,rfvars], type='prob')[,2]
	saveRDS(mapValidCells, file.path('dat', 'eval', paste0(spName, '_evalCells.rds')))
	
	# compute tss/auc and return
	with(mapValidCells, rbind(
		Find.Optim.Stat(Stat="TSS", Fit=(1000*fit.rf), Obs=obs),
		Find.Optim.Stat(Stat="ROC", Fit=(1000*fit.rf), Obs=obs)))
}, mc.cores=numCores)
names(eval.rf) = speciesList
saveRDS(eval.rf, 'res/eval/sdm_eval.rds')