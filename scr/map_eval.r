library(reshape2)
library(foreach)
library(doParallel)
library(raster)
library(rgdal)
library(randomForest)
library(biomod2)
registerDoParallel(cores=4)
source('scr/stm_functions.r')

pp = readRDS('dat/pp_grid_tall.rds')
tempPlotsAll = readRDS("dat/tempPlot_presence.rds")

# reshape to be wide by variable
ras.wide = dcast(pp, year_measured + lon + lat ~ biovar, value.var = 'val')
ras.wide = ras.wide[complete.cases(ras.wide),]

ras.wide[,'mean_diurnal_range'] = ras.wide[,'mean_diurnal_range'] / 10
ras.wide[,'mean_temp_wettest_quarter'] = ras.wide[,'mean_temp_wettest_quarter'] / 10

# scale variables
climScale = readRDS("dat/climate_scaling.rds")
ras.wide.scale = ras.wide
for(v in names(climScale$center))
	ras.wide.scale[,v] = scale(ras.wide.scale[,v], center = climScale$center[v], scale=climScale$scale[v])


# make a list by year
rasList = foreach(year=unique(ras.wide.scale$year_measured), .inorder=TRUE, 
		.final=function(x) {
			names(x) = unique(ras.wide.scale$year_measured)
			x}) %dopar%
{
	dfsub = ras.wide.scale[ras.wide.scale$year_measured == year,2:ncol(ras.wide.scale)]
	coordinates(dfsub) = 1:2
	gridded(dfsub) = TRUE
	brick(dfsub)
}


spName = '18032-ABI-BAL'



# read in calibration data and find points that were only sampled once
# (and thus not used to fit the STM)
trans = readRDS(file.path('dat', 'transition', paste0(spName, '_transitions.rds')))
pres = readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))
plotLocs = readRDS('dat/raw/plotLocations.rds')
singlePlots = pres[which(!(pres$plot_id %in% trans$plot)),]
singlePlots = merge(singlePlots, plotLocs)


# combine with temporary plots
tempPlots = tempPlotsAll[,c('plot_id', 'year_measured', spName, 'latitude', 'longitude')]
names(tempPlots) = names(singlePlots)
allPlots = rbind(singlePlots, tempPlots)
allPlots = allPlots[complete.cases(allPlots),]


# for each year, compute obs and converte to data frame with complete cases only
mapValidCells = do.call(rbind, lapply(unique(allPlots$year), function(year) {
	pl = allPlots[allPlots$year == year,c(spName, 'lon', 'lat')]
	coordinates(pl) = c('lon', 'lat')
	rb = rasList[[as.character(year)]]
	rb$obs = rasterize(pl, rb, 
			fun=function(x, ...) as.integer(sum(x, ...) > 0))[[2]]
	res = getValues(rb)
	res[complete.cases(res),]
}))


# get SDM fit
rfvars = c('annual_mean_temp', 'mean_diurnal_range', 'mean_temp_wettest_quarter', 
		'pp_seasonality', 'pp_warmest_quarter', 'tot_annual_pp')
RFMod = readRDS(file.path('res', 'sdm', paste0(spName, '_rf_sdm.rds')))
mapValidCells$fit.rf = predict(RFMod, newdata = mapValidCells[,rfvars], type='prob')[,2]
tss.rf = with(mapValidCells, Find.Optim.Stat(Stat="TSS", Fit=fit.rf, Obs=obs))

modName = '0'
sampleSize = NA

# get STM fit
posterior = readRDS(file.path('res', 'posterior', paste0(spName, '_posterior.rds')))[[modName]]
if(is.null(sampleSize) || is.na(sampleSize)) {
	samples = posterior
} else
	samples = posterior[sample(nrow(posterior), sampleSize),]


env1 = mapValidCells$annual_mean_temp
env2 = mapValidCells$tot_annual_pp
## mapValidCells$fit.stm = foreach(pars = iter(samples, by='row'), .combine=cbind, 
## 		.final=function(x) rowSums(x) / nrow(samples)) %dopar% {
## 	if(length(pars) == 2)
## 	{
## 		cp = pars[1]
## 		ep = pars[2]
## 	} else {
## 		cp = pars[1:5]
## 		ep = pars[6:10]
## 	}
## 	C = predict.stm_point(cp, env1, env2)
## 	E = predict.stm_point(ep, env1, env2)
## 	as.integer(C > E)
## }
# this is a single summary tss value
## tss.stm = with(mapValidCells, Find.Optim.Stat(Stat="TSS", Fit=fit.stm, Obs=obs))

system.time(tss.posterior <- foreach(pars = iter(samples, by='row'), .combine=c) %dopar% {
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
	fit = as.integer(C > E)
	Find.Optim.Stat(Stat="TSS", Fit=fit, Obs=mapValidCells$obs)[1]
}
)