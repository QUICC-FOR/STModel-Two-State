#!/usr/bin/env Rscript

# prepares temp plot data and evaluates the SDM/STM for all species specified
# gets posterior tss/roc at the map scale
# same as script C, except this script computes ROC for posterior replicates, rather than for the mean of the posterior
# yields a posterior distribution of ROC values
# script c MUST be completed for this to run successfully
# produces:
#    res/eval/stm_eval_posterior.rds

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



permPlotRasterFile = 'dat/clim/permPlot_annual_clim_raster.rds'
permPlot.annualClimRaster = readRDS(permPlotRasterFile)

tempPlotRasterFile = 'dat/clim/tempPlot_annual_clim_raster.rds'
tempPlot.annualClimRaster = readRDS(tempPlotRasterFile)

tempPlotsAll = readRDS("dat/tempPlot_presence.rds")

#### now evaluate the STM

modName = '0'
for(spName in speciesList) {
	cat("Starting species", spName, "... ")
	# get STM fit
	samples = readRDS(file.path('res', 'posterior', paste(spName, modName, 'samples.rds', sep='_')))
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]
	eval.rf = readRDS(file.path('res', 'eval', paste0(spName, '_sdm_eval.rds')))
	mapValidCells = eval.rf[['validation.data']]
	env1 = mapValidCells$annual_mean_temp
	env2 = mapValidCells$tot_annual_pp
	
	roc.reps = foreach(pars = iter(samples, by='row'), .packages='biomod2', 
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
		stm.fit = as.integer(C > E)
		roc = Find.Optim.Stat(Stat="ROC", Fit=(1000*stm.fit), Obs=mapValidCells$obs)
		c('roc'=roc[1], roc.thresh=roc[2])
	}
	saveRDS(roc.reps, file.path('res', 'eval', paste0(spName, '_stm_eval_posterior.rds')))
	cat("finished\n")
}
