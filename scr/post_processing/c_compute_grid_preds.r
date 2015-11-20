#!/usr/bin/env Rscript

## depends:
##    res/posterior/*

## makes:
##    res/posteriorGrid/*

library(coda)
speciesList = readRDS('dat/speciesList.rds')
models = c('0', 'i0', 'g', 'ig')
speciesInfo = read.csv('dat/speciesInfo.csv')
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

if(length(speciesList) == 0) stop("Error: no species specified")
cat("Will process posteriors for:\n")
for(sp in speciesList) cat("  ", sp, "\n")
climScale = readRDS('dat/climate_scaling.rds')
suppressWarnings(dir.create(file.path('res', 'posteriorGrid'), recursive=TRUE))

climGrid = readRDS(file.path('dat', 'climateGrid_scaled.rds'))
gr_env1 = climGrid$annual_mean_temp
gr_env2 = climGrid$tot_annual_pp

# the posterior columns containing extinction and colonization parameters
e_pars = function(pars)	if(length(pars) > 2) pars[6:10] else pars[2]
c_pars = function(pars)	if(length(pars) > 2) pars[1:5] else pars[1]


for(spName in speciesList)
{

	posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
	posterior = readRDS(posteriorFname)
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'stm_calib.rds')))

	for(mod in models)
	{

		curPosterior = posterior[[mod]]

		grPredict_e = grPredict_c = matrix(NA, nrow = nrow(curPosterior), ncol=length(gr_env1))
		for(i in 1:nrow(curPosterior))
		{
			pars = curPosterior[i,]
			grPredict_e[i,] = predict.stm_point(e_pars(pars), gr_env1, gr_env2)
			grPredict_c[i,] = predict.stm_point(c_pars(pars), gr_env1, gr_env2)
		}		
		grLambda = grPredict_c - grPredict_e
		grPres = (grLambda > 0) * 1

		saveRDS(list(coordinates = climGrid[,c('lon', 'lat')], E = grPredict_e, 
				C = grPredict_c, lambda = grLambda, stmPres = grPres), 
				file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
	}
}
