#!/usr/bin/env Rscript

## computes the status of the species in every cell on the map
## this is done as a posterior distribution, so takes a bit of time

## depends:
##    res/posterior/$species_$model_samples.rds
##    res/eval/$species_stm_eval.rds
##    res/eval/$species_sdm_eval.rds
##    res/sdm/$species_sdm_projection.rds
##    dat/clim/climateGrid_scaled.rds
##    dat/stm_calib/$species_stm_calib.rds



## makes:
##    
##    

library(coda)
library(foreach)
library(doParallel)
library(abind)
library(biomod2)
library(sp)
library(raster)
library(rgdal)

registerDoParallel(cores=detectCores())

speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
source('scr/stm_functions.r')
mapProj = readRDS("dat/map_projections.rds")

sampleSize = 10
model = '0'

# options are TSS or ROC
evalStat = 'ROC'

suppressWarnings(dir.create(file.path('res', 'rangemaps'), recursive=TRUE))
suppressWarnings(dir.create(file.path('res', 'areas'), recursive=TRUE))


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


# get projected climate map
climGrid = readRDS(file.path('dat', 'clim', 'climateGrid_scaled.rds'))
temperature = readRDS(file.path('dat', 'clim', 'climateGrid_unscaled.rds'))[,c('x','y', 'annual_mean_temp')]
gr_env1 = climGrid$annual_mean_temp
gr_env2 = climGrid$tot_annual_pp

# the posterior columns containing extinction and colonization parameters
e_pars = function(pars)	if(length(pars) > 2) pars[6:10] else pars[2]
c_pars = function(pars)	if(length(pars) > 2) pars[1:5] else pars[1]

rde.comp = function(sdm, stm)
{
	# presence is 1, colonization is 2, extinction is 3, absence is NA
	rde = (1 * (stm & sdm)) + 
			(2 * (stm & !sdm)) + 
			(3 * (!stm & sdm))
	# change absences to NA and subtract 1 from all categories
	rde[rde == 0] = NA 
	rde - 1
}

lam = function(T, P, pars)
{
	Cl = pars[1] + pars[2]*T + pars[3]*P + pars[4]*T^2 + pars[5]*P^2
	El = pars[6] + pars[7]*T + pars[8]*P + pars[9]*T^2 + pars[10]*P^2
	plogis(Cl) - plogis(El)
}

for(spName in speciesList)
{
	cat("Starting species", spName, "\n")
	posteriorFname = file.path('res', 'posterior', paste(spName, model, 'samples.rds', sep='_'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=sampleSize),]
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'_stm_calib.rds')))

	# get sdm threshold
	rfGrid = readRDS(file.path('res','sdm',paste0(spName, '_sdm_projection.rds')))
	rfEval = readRDS(paste0('res/eval/', spName, '_sdm_eval.rds'))
	sdmThreshold = rfEval[['validation.results']][toupper(evalStat), 'cutoff']/1000
	# both maple thresholds were a bit high and produced smaller than expected present ranges
	if(spName == '28728-ACE-RUB' | spName == '28731-ACE-SAC') sdmThreshold = 0.2

	# picea mariana and pinus banksiana threshold were a bit low
	if(spName == '183302-PIC-MAR') sdmThreshold = 0.15
	if(spName =='183319-PIN-BAN') sdmThreshold = 0.1
	rfGrid = readRDS(file.path('res','sdm',paste0(spName, '_sdm_projection.rds')))

	spGrid = within(rfGrid, sdm.pres <- as.integer(sdm >= sdmThreshold))

	if(spName == '183302-PIC-MAR' | spName == '183319-PIN-BAN')
		warmLims = readRDS(paste0("res/mask/", spName, "_temp_mask.rds"))
	
	stmGrid = foreach(pars = iter(samples, by='row'), i=icount(), .combine=function(...) abind(..., along=3),
			.multicombine=TRUE) %dopar%
	{
		## posterior computations
		rval = cbind(E=predict.stm_point(e_pars(pars), gr_env1, gr_env2), 
			C=predict.stm_point(c_pars(pars), gr_env1, gr_env2))

		# correct for model spec errors for 2 spp that lack northern range limits
		if(spName == '183302-PIC-MAR' | spName == '183319-PIN-BAN')
			rval = cbind(rval, mask = gr_env1 < warmLims[,i])		
		rval
	}
	lamGrid = stmGrid[,'C',] - stmGrid[,'E',]

	if(spName == '183302-PIC-MAR')
	{
		presGrid = (lamGrid > 0 & (is.na(stmGrid[,'mask',]) | stmGrid[,'mask',] == 1)) * 1
	} else if(spName == '183319-PIN-BAN')
	{
		# pinban still gets weird results with the mask (because of precip) so we use the SDM to help a bit
		presGrid = (lamGrid > 0 & (is.na(stmGrid[,'mask',]) | stmGrid[,'mask',] == 1) & 
				climGrid$y >= min(climGrid$y[spGrid$sdm.pres == 1])) * 1	
	} else {
		presGrid = (lamGrid > 0) * 1
	}
	cat("     ", spName, ": finished posterior*grid computation\n")
	
	##################
	#################
	# big objects created here
	#	stmGrid
	#	lamGrid
	#	presGrid


	# start of former script D
	# compute the range disequilibrium state for each cell in the grid
	# uses SDM and STM thresholds from the model evaluation to determine pres-abs cutoffs
	cat("     ", spName, ": computing range map\n")

	# get optimal thresholds for AUC/TSS using evaluation dataset
	stmEval = readRDS(paste0('res/eval/', spName, '_stm_eval.rds'))
	stmThreshold = stmEval[paste(tolower(evalStat), 'thresh', sep='.')]/1000
	

	# figure out presence in stm
	spGrid$stm = rowSums(presGrid, na.rm=TRUE)/ncol(presGrid)
	spGrid$stm.pres = as.integer(spGrid$stm >= stmThreshold)
	
	# save stats from the big posterior map
	spGrid$col = rowMeans(stmGrid[,'C',])
	spGrid$col.sd = apply(stmGrid[,'C',], 1, sd)
	spGrid$ext = rowMeans(stmGrid[,'E',])
	spGrid$ext.sd = apply(stmGrid[,'E',], 1, sd)
	spGrid$lambda = rowMeans(lamGrid)
	spGrid$lambda.sd = apply(lamGrid, 1, sd)

	## work out range disequilibrium in posterior
	spGrid = within(spGrid,
	{
		rde.present <- rowSums(presGrid  & spGrid$sdm.pres, na.rm=TRUE)/ncol(presGrid)
		rde.expand <- rowSums(presGrid & (!spGrid$sdm.pres), na.rm=TRUE)/ncol(presGrid)
		rde.contract <- rowSums((!presGrid) & spGrid$sdm.pres, na.rm=TRUE)/ncol(presGrid)
	})
	
	# "consensus" range disequilibrium
	spGrid$rde = rde.comp(spGrid$sdm.pres, spGrid$stm.pres)	
	saveRDS(spGrid, file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))


	
	# start of former script E
	# determines the potential change in every species' range
	cat("     ", spName, ": calculating areas\n")

	area.sp = foreach(stm.pres = iter(presGrid, by='column'), .combine=rbind, 
			.final = function(x) {
				x = as.data.frame(x)
				colnames(x) = c('presence', 'expansion', 'extinction')
				x}) %dopar%
	{
		rde = cbind(rde.comp(spGrid$sdm.pres, stm.pres), spGrid[,c('x', 'y')])
		coordinates(rde) = c('x', 'y')
		gridded(rde) = TRUE
		proj4string(rde) = mapProj$projected
		rde = raster(rde)
		rde.freq = freq(rde)
		if(!(1 %in% rde.freq[,1])) rde.freq = rbind(rde.freq, c(1,0))
		if(!(2 %in% rde.freq[,1])) rde.freq = rbind(rde.freq, c(2,0))
		if(!(3 %in% rde.freq[,1])) rde.freq = rbind(rde.freq, c(3,0))
		rde.freq = rde.freq[order(rde.freq[,1]),]
		
		# the raster is in meters, we present results in 1000s of km2
		rde.freq[1:3,2] * prod(res(rde)/1000)/1000
	}
	saveRDS(area.sp, file.path('res','areas',paste0(spName, '_', model, '_areas.rds')))
}
