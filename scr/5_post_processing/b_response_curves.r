#!/usr/bin/env Rscript
# compute the posterior distribution of the response curves

## depends:
##    res/posterior/*

## makes:
##    res/resp_curve/*

library(coda)

# number of points in horizontal dimension of response curve
# number of posterior rows to use
rcRes = 1000
posterior.n = 1000


speciesList = readRDS('dat/speciesList.rds')
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
climScale = readRDS('dat/clim/climate_scaling.rds')
climDat = readRDS('dat/clim/plotClimate_scaled.rds')
suppressWarnings(dir.create(file.path('res', 'resp_curve'), recursive=TRUE))

# the posterior columns containing extinction and colonization parameters
e_pars = function(pars)	if(length(pars) > 2) pars[6:10] else pars[2]
c_pars = function(pars)	if(length(pars) > 2) pars[1:5] else pars[1]

# set ranges and non-constant x-values for the response curves
trange = scale(c(-5, 25), center = climScale$center['annual_mean_temp'], 
			scale = climScale$scale['annual_mean_temp'])
prange = scale(c(500, 2000), center = climScale$center['tot_annual_pp'], 
			scale = climScale$scale['tot_annual_pp'])
rc.env1 = seq(trange[1], trange[2], length.out=rcRes)
rc.env2 = seq(prange[1], prange[2], length.out=rcRes)

for(spName in speciesList)
{
	cat("starting", spName, "\n")
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_0_samples.rds'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'_stm_calib.rds')))

	# get the median precip and temp values where the species is present
	# for setting constants for response curves
	presDat = readRDS(file.path('dat', 'presence', paste0(spName,'_presence.rds')))
	presDat = merge(presDat, climDat)
	rc.env1.c = info$rc_tval
	rc.env2.c = info$rc_pval
	if(is.na(rc.env1.c))
		rc.env1.c = mean(presDat[presDat[,spName] == 1,'annual_mean_temp'])
	if(is.na(rc.env2.c))
		rc.env2.c = mean(presDat[presDat[,spName] == 1,'tot_annual_pp'])
	
	rcPredict1_e = rcPredict1_c = rcPredict2_e = rcPredict2_c = 
				matrix(nrow = nrow(samples), ncol=rcRes)
	for(i in 1:nrow(samples))
	{
		pars = samples[i,]
		rcPredict1_e[i,] = predict.stm_point(e_pars(pars), rc.env1, rc.env2.c)
		rcPredict1_c[i,] = predict.stm_point(c_pars(pars), rc.env1, rc.env2.c)
		rcPredict2_e[i,] = predict.stm_point(e_pars(pars), rc.env1.c, rc.env2)
		rcPredict2_c[i,] = predict.stm_point(c_pars(pars), rc.env1.c, rc.env2)
	}
		
	respCurve = data.frame(
		temp = (rc.env1 * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp'],
		col.temp = colMeans(rcPredict1_c),
		col.temp.upper = apply(rcPredict1_c,2,quantile,0.95),
		col.temp.lower = apply(rcPredict1_c,2,quantile,0.05),
		ext.temp = colMeans(rcPredict1_e),
		ext.temp.upper = apply(rcPredict1_e,2,quantile,0.95),
		ext.temp.lower = apply(rcPredict1_e,2,quantile,0.05),
		precip = (rc.env2 * climScale$scale['tot_annual_pp']) + climScale$center['tot_annual_pp'],
		col.precip = colMeans(rcPredict2_c),
		col.precip.upper = apply(rcPredict2_c,2,quantile,0.95),
		col.precip.lower = apply(rcPredict2_c,2,quantile,0.05),
		ext.precip = colMeans(rcPredict2_e),
		ext.precip.upper = apply(rcPredict2_e,2,quantile,0.95),
		ext.precip.lower = apply(rcPredict2_e,2,quantile,0.05)
	)
	saveRDS(respCurve, file.path('res','resp_curve',paste0(spName,'_0_respCurve.rds')))

}


