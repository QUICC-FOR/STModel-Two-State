#!/usr/bin/env Rscript

# Pinus banksiana and Picea mariana have response curves that aren't informed enough
# to understand that there is a missing range limit in the north
# these curves bend creating a pseudo second southern re-occurrence limit, producing a weird
# map with presence, then a hole, then new presence (moving north to south)
# to fisk this, we need to find the real southern range limit via looking for the 
# root of the response curve

library(doParallel)
library(foreach)
library(coda)
registerDoParallel(cores=detectCores())

# get projected climate map
climGrid = readRDS(file.path('dat', 'clim', 'climateGrid_scaled.rds'))
temperature = readRDS(file.path('dat', 'clim', 'climateGrid_unscaled.rds'))[,c('x','y', 'annual_mean_temp')]
gr_env1 = climGrid$annual_mean_temp
gr_env2 = climGrid$tot_annual_pp

lam = function(T, P, pars)
{
	Cl = pars[1] + pars[2]*T + pars[3]*P + pars[4]*T^2 + pars[5]*P^2
	El = pars[6] + pars[7]*T + pars[8]*P + pars[9]*T^2 + pars[10]*P^2
	plogis(Cl) - plogis(El)
}


library(rootSolve)
tmin = min(gr_env1)
tmax = max(gr_env1)
getroot = function(p, pars)
{
	root = uniroot.all(lam, P=p, pars=pars, lower=tmin, upper=tmax)
	if(length(root) > 1)
		root = root[which(lam(root-0.05, p, pars) > 0 & lam(root + 0.05, p, pars) < 0)]
	if(length(root) == 0)
		root = NA
	root
}

speciesList = c("183302-PIC-MAR", "183319-PIN-BAN")
for(spName in speciesList)
{
	posteriorFname = file.path('res', 'posterior', paste(spName, '0', 'samples.rds', sep='_'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=1000),]
	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))

	# remove unscaled and add scaled temperature
	spGrid = spGrid[,-which(colnames(spGrid) == 'annual_mean_temp')]
	spGrid = merge(spGrid, climGrid[,c('x', 'y', 'annual_mean_temp', 'tot_annual_pp')])
	
	
	tooWarm = foreach(pars = iter(samples, by='row'), .combine=cbind) %dopar%
	{
		# find the southern range limit
		# we do this by solving for the root, starting from the coldest temperature
		sapply(gr_env2, getroot, pars=pars)
	}	
	saveRDS(tooWarm, paste0('res/', spName, 'temp_mask.rds'))
			
}			
