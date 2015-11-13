# process temp plot climate
climScale = readRDS("dat/climate_scaling.rds")
tp_clim = read.csv("dat/raw/temp_plots/tp_climData_reshaped.csv", dec=',', sep=';')
colsToKeep = c('plot_id', 'year_measured')
tp_clim = cbind(tp_clim[,colsToKeep], scale(tp_clim[,names(climScale$center)], 
		center=climScale$center, scale=climScale$scale))
tp_clim = tp_clim[complete.cases(tp_clim),]


tempPlots = read.csv("dat/raw/temp_plots/tp_plotinfoData.csv", dec='.', sep=';')

tp_species = read.csv("dat/raw/temp_plots/tp_treeData_allSpecies.csv", dec='.', sep=';')

treeDat = merge(tp_species, tempPlots, all.x=TRUE)

# drop unneeded stuff
treeDat = treeDat[-(which(treeDat$id_spe == "")),-7]
sampleDat = dcast(treeDat, plot_id + year_measured + latitude + longitude ~ id_spe, fill = 0, 
		value.var = "basal_area", fun.aggregate = function(x) as.integer(sum(x) > 0))
stateData = merge(sampleDat, tp_clim, all = 'T', by=c("plot_id", "year_measured"))
## stateData = merge(stateData, tempPlots, all = 'T', by=c("plot_id", "year_measured"))
stateData = stateData[complete.cases(stateData),]
saveRDS(stateData, "dat/tempPlot_presence.rds")