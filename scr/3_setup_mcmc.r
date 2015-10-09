#!/usr/bin/Rscript

mcmcProportion = 0.75
stmMaxInterval = 15
env1 = "annual_mean_temp"
env2 = "tot_annual_pp"
source("stm_functions.r")
dir.create(file.path('dat', 'stm_calib'), recursive=TRUE)
dir.create(file.path('dat', 'stm_valid'), recursive=TRUE)
dir.create(file.path('dat', 'mcmc'), recursive=TRUE)

plotLocations = readRDS('dat/raw/plotLocations.rds')
transitionClimate = readRDS('dat/transitionClimate_scaled.rds')
transitionClimate = merge(transitionClimate, plotLocations, by.x=1, by.y=1)


subset_transitions = function(dat, frac=0.5, max.interval=15, mask.tol = c(10,10))
{
	# subset the transition data
	# first, apply the mask and drop all intervals greater than 15 years
	# then take half (or whatever fraction) of the remaining observed transitions
	# along with the same fraction of non-transitions
	intervals = dat$year2 - dat$year1
	dat.mask = stm_mask(dat[,c('lon', 'lat')], as.integer(dat$state1 | dat$state2), 
			dat[,c('lon', 'lat')], mask.tol)
	indices = which(intervals <= max.interval & dat.mask == 1)
	datMasked = dat[indices,]
	
	# find the indices for each of the 4 transition types
	trIndices = with(datMasked, list(
		col = which(state1 == 0 & state2 == 1),
		ext = which(state1 == 1 & state2 == 0),
		pre = which(state1 == 1 & state2 == 1),
		abs = which(state1 == 0 & state2 == 0)))
	
	# select an equal fraction of each type
	sel = sapply(trIndices, function(x) sample(x, ceiling(frac*length(x))))
	
	list(selected = datMasked[sel,], unselected = datMasked[-sel,])
}



# Loop across species

for(spName in speciesList)
{
	# read in transition data and split into calibration and validation sets
	# be sure to filter by the STM Mask BEFORE splitting - data outside the mask is dropped, not reserved
	transitionAll = readRDS(file.path("dat", "transition", paste0(spName, "_transitions.rds")))
	sdm = readRDS(file.path('res', 'sdm', paste0(spName, "_rf_sdm.rds")))
	prevalence = predict(sdm, newdata = transitionClimate, type='prob')[,2]
	climCols = c('plot', 'lon', 'lat', 'year1', 'year2', env1, env2)
	spTrClim = transitionClimate[,climCols]
	transitionAll = merge(transitionAll, spTrClim)
	transitionAll$prevalence = prevalence

	transitionSubsets = subset_transitions(transitionAll, mcmcProportion, stmMaxInterval, stmMaskTolerance)
	saveRDS(transitionSubsets$selected, file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
	saveRDS(transitionSubsets$unselected, file.path('dat', 'stm_valid', paste0(spName, 'stm_valid.rds')))

	# prep the mcmc transition file
	transitions = transitionSubsets$selected
	mcmcDataFile = file.path('dat', 'mcmc', paste0(spName, '_mcmcDat.txt'))
	mcmcData = data.frame(
		initial = transitions$state1,
		final = transitions$state2,
		env1 = transitions[,env1],
		env2 = transitions[,env2],
		interval = transitions$year2-transitions$year1,
		prevalence1 = transitions$prevalence)
	write.csv(mcmcData, mcmcDataFile, row.names=FALSE)


	# prep the inits
	parNames = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
	isConstant = rep(0,length(parNames))
	isConstant[parNames %in% c('g5','g6','e5','e6')] = 1
	mcmcInitFile = file.path('dat', 'mcmc', paste0(spName, '_mcmcInit.txt'))
	mcmcInits = data.frame(
		name = parNames,
		initialValue = 0,
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = isConstant)
	mcmcInits$priorSD[parNames %in% c('g0', 'e0')] = 10
	write.csv(mcmcInits, mcmcInitFile, row.names=FALSE)
}