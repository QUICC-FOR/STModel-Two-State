#!/usr/bin/env Rscript

# depends:
#    res/posterior/*
#    res/eval/ras/*

# produces:
#    res/eval/SPNAME_MOD__rdeAreas.rds

library(foreach)
library(doParallel)
library(coda)
library(raster)
library(sp)
source('scr/stm_functions.r')


## clArgs <- commandArgs(trailingOnly = TRUE)
## 
## ## these will be read in from the command line
## spName = clArgs[1]
## modName = if(length(clArgs) > 1) clArgs[2] else '0'
## numCores = if(length(clArgs) > 2) clArgs[3] else detectCores()
## sampleSize = if(length(clArgs) > 3) clArgs[4] else NA
## 

speciesList = readRDS('dat/speciesList.rds')
modName = '0'
numCores = 7
sampleSize = 1000
registerDoParallel(cores=numCores)

for(spName in speciesList)
{

cat('starting', spName, '\n')

	gr_env1 = readRDS(file.path('res', 'eval', 'ras', paste0(spName, '_env1.rds'))) 
	gr_env2 = readRDS(file.path('res', 'eval', 'ras', paste0(spName, '_env1.rds'))) 
	sdmPres = readRDS(file.path('res', 'eval', 'ras', paste0(spName, '_sdmPres.rds'))) 
	stmProb = readRDS(file.path('res', 'eval', 'ras', paste0(spName, '_stmProb.rds'))) 


	rdeCats = list('absent' = 0, 'present' = 1, 'expand' = 2, 'contract' = 3)

	# grab the posterior samples
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
	posterior = readRDS(posteriorFname)[[modName]]

	if(is.null(sampleSize) || is.na(sampleSize)) {
		samples = posterior
	} else
		samples = posterior[sample(nrow(posterior), sampleSize),]


	# foreach posterior replicate
	rdePosterior = foreach(curPars = iter(samples, by='row'), .combine=rbind, 
					.packages = 'raster', .final=function(x) {
		colnames(x) = names(rdeCats)[2:4]
		x}) %dopar%
	{
		# 1. compute C & E (as rasters! -- this saves having to reproject at each rep)
		if(length(curPars) == 2)
		{
			cPars = curPars[1]
			ePars = curPars[2]
		} else {
			cPars = curPars[1:5]
			ePars = curPars[6:10]
		}
		C = predict.stm_point(cPars, gr_env1, gr_env2)
		E = predict.stm_point(ePars, gr_env1, gr_env2)

		# 2. compute STM presence
## note - this was a bit of a challenge
## the first way, many areas end up with C>E but with VERY small rates
## so we probably overestimate the area that's colonized
## the second way might be a bit too agressive the other direction, but the values
## more closely match the map
## 		stmPres = (C > E)
		stmPres = (C > E) & stmProb > 0.1
	
		# 3. compute RDE
		rde = (rdeCats[['present']] * (stmPres & sdmPres)) + 
			(rdeCats[['expand']] * (stmPres & !sdmPres)) + 
			(rdeCats[['contract']] * (!stmPres & sdmPres))
		# change absences to NA and subtract 1 from all cats
		rde[rde==0] = NA
		rde = rde - 1
	
		# 5. compute zonal area
		cellArea_km2 = prod(res(rde)/1000)
		freqs = freq(rde)
		(freqs[,2] * cellArea_km2)[1:3]
	}

	outFile = file.path('res','eval',paste0(spName, '_', modName, '_rdeAreas.rds'))
	saveRDS(rdePosterior, outFile)
}