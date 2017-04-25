#!/usr/bin/env Rscript
library("INLA")

source("scr/inla/inla_functions.r")

spAll <- readRDS("dat/speciesList.rds")
spList <- commandArgs()
spList <- spList[spList %in% spAll]
if(length(spList) == 0) spList <- spAll

# min and max intervals to use
minInterval <- 5
maxInterval <- 5

## spName <- "28731-ACE-SAC"


climGrid <- readRDS("dat/clim/climateGrid_scaled.rds")
coords.pr <- as.matrix(climGrid[,c('x', 'y')])



for(spName in spList)
{
	resDir <- file.path('res', 'inla', spName)
	dir.create(resDir, showWarnings = FALSE)

	# calibration data, with extinctions pulled out
	dat <- readRDS(paste0("dat/stm_calib/", spName, "_stm_calib.rds"))
	dat$interval <- dat$year2 - dat$year1
	dat <- dat[dat$interval >= minInterval & dat$interval <= maxInterval & dat$state1 == 1,]
	dat$ext <- with(dat, as.integer(state1 & !(state2)))

	message(Sys.time(), " starting extinction model for ", spName)
	emod <- ext_fit(dat)
	saveRDS(emod, file.path(resDir, "inla_models_ext.rds"))

	
	message("    ", Sys.time(), " sampling from model")
	samps <- get_inla_samples(emod, coords.pr)
	saveRDS(samps, file.path(resDir, "inla_samples.rds"))
	
	
	message("    ", Sys.time(), " projecting model")
	mcmc <- read.csv(file.path('res', '/mcmc/', spName, '0', '/posterior.csv'))
	mcmc <- mcmc[seq(1, nrow(mcmc), length.out=1000),8:12]  ## extinction parameters only
	lp <- inla_lp(samps$samples, mcmc, samps$A.mat, climGrid$annual_mean_temp, climGrid$tot_annual_pp)
	lpGrid <- grid_summary(lp, coords.pr)
	saveRDS(lpGrid, file.path(resDir, 'lin_pred_grid.rds'))
	
	# apply reverse of link function
	exRates <- lapply(lp, plogis)
	exGrid <- grid_summary(exRates, coords.pr)
	saveRDS(exGrid, file.path(resDir, 'ex_rate_grid.rds'))
}

