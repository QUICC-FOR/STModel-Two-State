library(foreach)
library(doParallel)
library(reshape2)
library(randomForest)
library(coda)
registerDoParallel(cores=4)
## setwd("/Users/mtalluto/Dropbox/work/projects/STModel-Two-State_git/")


## notes to self
### this isn't really going to work the way it's done here. The temporary plot data
# are only for quebec. this means there isn't much variance in climate, and there is lots
# of spatial duplication and poor resolution in the STM predictions (mostly just 0 or 1)
# need to read in the original dataset and get any calibration data points that are NOT in the
# transition dataset (i.e., plots that were only sampled once)
# then need to intersect them with the climate rasters, so that the observed values are 1 or 0 for
# each climate cell (within a given year)
# the output will be a set of rasters for each year, with values of 1 (at least one observed presence), 
# 0 (no presences), or NA (no observations for that cell for that year)



source("scr/stm_functions.r")
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
sampleDat = dcast(treeDat, plot_id + year_measured + latitude + longitude ~ id_spe, fill = 0, 
		value.var = "basal_area", fun.aggregate = function(x) as.integer(sum(x) > 0))
stateData = merge(sampleDat, tp_clim, all = 'T', by=c("plot_id", "year_measured"))
## stateData = merge(stateData, tempPlots, all = 'T', by=c("plot_id", "year_measured"))
stateData = stateData[complete.cases(stateData),]

spName = '18032-ABI-BAL'
modName = '0'

# to do:
# √1. for grid cells with multiple plots, select one at random
# √2. read in the original random forest model (the one used to fit the STM) and project to temp plots
# 3. read STM pars and project to temporary plot grid (one per posterior rep)
# 4. compute TSS for RF and for STM for both datasets

# 1. for grid cells with multiple plots, select one at random
# get unique IDs based on climate values and select one plot at random from each set of identical climate values
climCols = c(colnames(tp_clim), spName, 'longitude', 'latitude')
climDat = stateData[,climCols]
# subset climDat to the projection range (+/- 10 degrees of the most extreme presences)
trans = readRDS(file.path("dat", "transition", paste0(spName, "_transitions.rds")))
plotLocations = readRDS('dat/raw/plotLocations.rds')
trans = merge(trans, plotLocations, all=TRUE, by.x = 'plot', by.y = 'plot_id')
trans = trans[complete.cases(trans),]
rows = stm_mask(new.coords = climDat[,c('longitude', 'latitude')], 
	pres = as.integer(trans$state1 | trans$state2), pres.coords = trans[,c('lon', 'lat')])
climDat = climDat[rows==1,]


rfvars = c('annual_mean_temp', 'mean_diurnal_range', 'mean_temp_wettest_quarter', 
		'pp_seasonality', 'pp_warmest_quarter', 'tot_annual_pp')
climDat.unique = unique(climDat[,rfvars])
climDat.unique$climID = 1:nrow(climDat.unique)
climDat = merge(climDat, climDat.unique)
climIDID = with(climDat, ave(annual_mean_temp, factor(climID), FUN=function(x) sample.int(length(x))))
climDat.subset = climDat[climIDID == 1,]

# 2. read in the original random forest model (the one used to fit the STM) and project to temp plots
RFMod = readRDS(file.path('res', 'sdm', paste0(spName, '_rf_sdm.rds')))
rfvars = c('annual_mean_temp', 'mean_diurnal_range', 'mean_temp_wettest_quarter', 
		'pp_seasonality', 'pp_warmest_quarter', 'tot_annual_pp')
climDat$rf.pred = predict(RFMod, newdata = climDat[,rfvars], type='prob')[,2]
climDat.subset$rf.pred = predict(RFMod, newdata = climDat.subset[,rfvars], type='prob')[,2]

# 3. read STM pars and project to temporary plot grid (one per posterior rep), compute TSS
## stmvars = c('annual_mean_temp', 'tot_annual_pp')
## cdat = climDat
cdat = climDat.subset
## 
STMPars = readRDS(file.path('res', 'posterior', paste0(spName, '_posterior.rds')))[[modName]]
STMPars = STMPars[sample(nrow(STMPars), 1000),]
obs = cdat[,spName]
numPres = foreach(pars = iter(STMPars, by='row'), .combine=cbind, .final=function(x) rowSums(x)) %dopar% {
## foreach(pars = t(STMPars), .combine=c, .packages='biomod2') %do% {
	if(length(pars) == 2) {
		C = predict.stm_point(pars[1], env1 = cdat[,'annual_mean_temp'], env2 = cdat[,'tot_annual_pp'])
		E = predict.stm_point(pars[2], env1 = cdat[,'annual_mean_temp'], env2 = cdat[,'tot_annual_pp'])
	} else {
		C = predict.stm_point(pars[1:5], env1 = cdat[,'annual_mean_temp'], env2 = cdat[,'tot_annual_pp'])
		E = predict.stm_point(pars[6:10], env1 = cdat[,'annual_mean_temp'], env2 = cdat[,'tot_annual_pp'])
	}
	as.integer(C > E)
}


colFun = function(x) rgb(colorRamp(c('black', 'red'))(x), maxColorValue=256)
plot(cdat$longitude, cdat$latitude, col=colFun(cdat$rf.pred), pch=16, cex=0.5)

quartz()
plot(cdat$longitude, cdat$latitude, col=colFun(numPres/nrow(STMPars)), pch=16, cex=0.5)

# 6. Read in original calibration P/A dataset
# 7. Aggregate PA to 1 or 0 for each climate grid cell for calibration and temp plot data
# 8. Build a new random forest on these aggregate data for the calibration set
# 9. Project RF and STM to the aggregated temporary plot dataset
# 10. Compute TSS for both models for the aggregated temp plot dataset