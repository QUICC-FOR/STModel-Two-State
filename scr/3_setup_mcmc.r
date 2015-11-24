#!/usr/bin/env Rscript

## prepare the data files for the MCMC
## sets up an initial run of the sampler with null inits; this should be run for just a 
## few samples to get an idea of inits to use

## dependencies:
##     dat/raw/plotLocations.rds
##     dat/transitionClimate_scaled.rds (script 1)
##     dat/transition/*
##     res/sdm/* (script 2)

## creates:
##     dat/stm_calib/*
##     dat/mcmc/*_mcmcDat.txt
##     dat/mcmc/*_mcmcInit.txt
##     scr/mcmc/*

library(randomForest)
## mcmcProportion = 1
stmMaxInterval = 15
env1 = "annual_mean_temp"
env2 = "tot_annual_pp"
source("scr/stm_functions.r")

# dir.create warns you if the directory already exists
suppressWarnings(
{
	dir.create(file.path('dat', 'stm_calib'), recursive=TRUE)
	dir.create(file.path('dat', 'mcmc'), recursive=TRUE)
	dir.create(file.path('res', 'mcmc'), recursive=TRUE)
	dir.create(file.path("log"), recursive=TRUE)
})
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

plotLocations = readRDS('dat/raw/plotLocations.rds')
transitionClimate = readRDS('dat/clim/transitionClimate_scaled.rds')
transitionClimate = merge(transitionClimate, plotLocations, by.x=1, by.y=1)


subset_transitions = function(dat, max.interval=15, mask.tol = c(10,10))
{
	# apply the mask and drop all intervals greater than 15 years
	intervals = dat$year2 - dat$year1
	dat.mask = stm_mask(dat[,c('lon', 'lat')], as.integer(dat$state1 | dat$state2), 
			dat[,c('lon', 'lat')], mask.tol)
	indices = which(intervals <= max.interval & dat.mask == 1)
	dat[indices,]
}



# Loop across species

for(spName in speciesList)
{
	# read in transition data and process into calibration set
	cat(paste("setting up species", spName, "\n"))
	transitionAll = readRDS(file.path("dat", "transition", paste0(spName, "_transitions.rds")))
	sdm = readRDS(file.path('res', 'sdm', paste0(spName, "_rf_sdm.rds")))
	prevalence = predict(sdm, newdata = transitionClimate, type='prob')[,2]
	climCols = c('plot', 'lon', 'lat', 'x', 'y', 'year1', 'year2', env1, env2)
	spTrClim = transitionClimate[,climCols]
	transitionAll = merge(transitionAll, spTrClim)
	transitionAll$prevalence = prevalence

	transitions = subset_transitions(transitionAll, stmMaxInterval, stmMaskTolerance)
	saveRDS(transitions, file.path('dat', 'stm_calib', paste0(spName, '_stm_calib.rds')))

	# prep the mcmc transition file
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
	mcmcIntInitFile = file.path('dat', 'mcmc', paste0(spName, '_mcmcInit_int.txt'))
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
	
	mcmcIntInits = mcmcInits
	mcmcIntInits$isConstant[!(parNames %in% c('g0', 'e0'))] = 1
	write.csv(mcmcIntInits, mcmcIntInitFile, row.names=FALSE)
}