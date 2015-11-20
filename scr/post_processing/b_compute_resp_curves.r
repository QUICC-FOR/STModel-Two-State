#!/usr/bin/env Rscript

## depends:
##    res/posterior/*

## makes:
##    res/resp_curve/*

library(coda)
library(doParallel)



numCores=detectCores() - 2
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
		
	corTest = as.numeric(arg[length(arg)])
	if(!is.na(corTest)) numCores = corTest
}
registerDoParallel(cores=numCores)

if(length(speciesList) == 0) stop("Error: no species specified")
cat("Will process posteriors for:\n")
for(sp in speciesList) cat("  ", sp, "\n")
climScale = readRDS('dat/climate_scaling.rds')
suppressWarnings(dir.create(file.path('res', 'resp_curve'), recursive=TRUE))




# the posterior columns containing extinction and colonization parameters
e_pars = function(pars)	if(length(pars) > 2) pars[6:10] else pars[2]
c_pars = function(pars)	if(length(pars) > 2) pars[1:5] else pars[1]

# number of points in horizontal dimension of response curve
rcRes = 1000

for(spName in speciesList)
{

	posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
	posterior = readRDS(posteriorFname)
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'stm_calib.rds')))

	for(mod in models)
	{

		curPosterior = posterior[[mod]]

		# set x values for the response curves; for constants (env1 is constant in the curve 
		# for env2, we lookup the values in info; if missing we use 0 (the mean)
		trange = scale(c(-5, 25), center = climScale$center['annual_mean_temp'], 
					scale = climScale$scale['annual_mean_temp'])
		prange = scale(range(calibDat$tot_annual_pp), 
					center = climScale$center['tot_annual_pp'], 
					scale = climScale$scale['tot_annual_pp'])
		rc_env1 = seq(trange[1], trange[2], length.out=rcRes)
		rc_env1_c = if(is.null(info$rc_tval) || is.na(info$rc_tval)) 
			{0} else {info$rc_tval}
		rc_env2 = seq(prange[1], prange[2], length.out=rcRes)
		rc_env2_c = if(is.null(info$rc_pval) || is.na(info$rc_pval))
			{0} else {info$rc_pval}
			
		rcPredict1_e = foreach(pars = iter(curPosterior, by='row'), .combine=rbind) %dopar%
			predict.stm_point(e_pars(pars), rc_env1, rc_env2_c)
		rcPredict1_c = foreach(pars = iter(curPosterior, by='row'), .combine=rbind) %dopar%
			predict.stm_point(c_pars(pars), rc_env1, rc_env2_c)
		rcPredict2_e = foreach(pars = iter(curPosterior, by='row'), .combine=rbind) %dopar%
			predict.stm_point(e_pars(pars), rc_env1_c, rc_env2)
		rcPredict2_c = foreach(pars = iter(curPosterior, by='row'), .combine=rbind) %dopar%
			predict.stm_point(c_pars(pars), rc_env1_c, rc_env2)

		respCurve = data.frame(
			temp = (rc_env1 * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp'],
			col.temp = colMeans(rcPredict1_c),
			col.temp.upper = apply(rcPredict1_c,2,quantile,0.95),
			col.temp.lower = apply(rcPredict1_c,2,quantile,0.05),
			ext.temp = colMeans(rcPredict1_e),
			ext.temp.upper = apply(rcPredict1_e,2,quantile,0.95),
			ext.temp.lower = apply(rcPredict1_e,2,quantile,0.05),
			precip = (rc_env2 * climScale$scale['tot_annual_pp']) + climScale$center['tot_annual_pp'],
			col.precip = colMeans(rcPredict2_c),
			col.precip.upper = apply(rcPredict2_c,2,quantile,0.95),
			col.precip.lower = apply(rcPredict2_c,2,quantile,0.05),
			ext.precip = colMeans(rcPredict2_e),
			ext.precip.upper = apply(rcPredict2_e,2,quantile,0.95),
			ext.precip.lower = apply(rcPredict2_e,2,quantile,0.05)
		)
		saveRDS(respCurve, file.path('res','resp_curve',paste0(spName,'_', mod, '_respCurve.rds')))

	}
}


