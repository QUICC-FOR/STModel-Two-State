#!/usr/bin/env Rscript

## reads all posterior results in and processes them into CODA objects; also makes some plots
## depends:
##    res/mcmc/*

## makes:
##    res/posterior/*
##    img/mcmc/*

library(coda)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
source('scr/stm_functions.r')
climScale = readRDS('dat/climate_scaling.rds')
suppressWarnings(dir.create(file.path('img', 'mcmc'), recursive=TRUE))
suppressWarnings(dir.create(file.path('res', 'posterior'), recursive=TRUE))
suppressWarnings(dir.create(file.path('res', 'resp_curve'), recursive=TRUE))
suppressWarnings(dir.create(file.path('res', 'maps'), recursive=TRUE))

climGrid = readRDS(file.path('dat', 'climateGrid_scaled.rds'))
gr_env1 = climGrid$annual_mean_temp
gr_env2 = climGrid$tot_annual_pp

sdmThreshold = 0.1


remove_constant_cols = function(x, tol=1e-10)
{
	vars = apply(x, 2, var)
	cols = which(vars > tol)
	x[,cols]
}

# the posterior columns containing extinction and colonization parameters
eCols = 6:10
cCols = 1:5

plot_mcmc_summary = function(PD, sp, density = FALSE, make.png = TRUE)
{
	dens = FALSE
	tr = TRUE
	w = 10
	h = 15
	dpi=600
	fontsize=11
	fname = file.path('img','mcmc', paste0(spName, '_traceplots.png'))
	
	if(density)
	{
		dens = TRUE
		tr = FALSE
		fname = file.path('img','mcmc', paste0(spName, '_densplots.png'))
	}
	if(make.png)
	{
		png(width=as.integer(dpi*w), height=as.integer(dpi*h), file=fname, 
				pointsize=fontsize, res=dpi)
	} else
	{
		quartz(width = w, height=h)
	}
	par(mfcol=c(10,4), mar=c(3,3,0,0), mgp=c(1.5, .5, 0), oma=c(0,0,2,2), cex=0.5, bty='n')
	for(mod in PD)
	{
		if(ncol(mod) == 2) plot(0,0,type='n', xaxt='n', yaxt='n', xlab='', ylab='')
		plot(mod, density=dens, trace=tr, ask=FALSE, auto.layout=FALSE, main='')
	}
	if(make.png) dev.off()
}

posteriors = list()
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	models = list.files(file.path('res','mcmc',spName))
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
	
	## check to see if posterior distributions for the species exist
	spPosterior = ifelse(file.exists(posteriorFname), 
			readRDS(posteriorFname), list())
	for(mod in models)
	{
		pos = read.csv(file.path('res','mcmc',spName,mod,'posterior.csv'))
		pos = remove_constant_cols(pos)
		pd = mcmc(pos, start=info$thin*info$burnin+1, thin=info$thin)
		spPosterior[[mod]] = pd
	}
	# drop any pesky null elements created by list() calls
	spPosterior	= spPosterior[-sapply(spPosterior, is.null)]
	saveRDS(spPosterior, posteriorFname)
	posteriors[[spName]] = spPosterior

	# set up traceplots and density plots
	plot_mcmc_summary(spPosterior, spName)
	plot_mcmc_summary(spPosterior, spName, density=TRUE)	
}

# the following will set up the DATA for plotting posterior distributions
# but it will not actually make the plots 
# plotting is fast, tweaking data is slow; this is a slow script

# number of points in horizontal dimension of response curve
rcRes = 1000
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'stm_calib.rds')))

	# set x values for the response curves; for constants (env1 is constant in the curve 
	# for env2, we lookup the values in info; if missing we use 0 (the mean)
	rc_env1 = with(calibDat, seq(min(annual_mean_temp), max(annual_mean_temp), length.out=rcRes))
	rc_env1_c = if(is.null(info$rc_tval) || is.na(info$rc_tval)) 
		{0} else {info$rc_tval}
	rc_env2 = with(calibDat, seq(min(tot_annual_pp), max(tot_annual_pp), length.out=rcRes))
	rc_env2_c = if(is.null(info$rc_pval) || is.na(info$rc_pval))
		{0} else {info$rc_pval}
	rcPredict1_e = rcPredict1_c = rcPredict2_e = rcPredict2_c = matrix(NA, nrow=nrow(spPosterior), ncol=rcRes)
	grPredict_e = grPredict_c = grLambda = grPres = matrix(NA, nrow=nrow(spPosterior), ncol=length(gr_env1))
	for(i in 1:nrow(spPosterior))
	{
		# response curve predictions
		rcPredict1_e[i,] = predict.stm_point(spPosterior[i,eCols], rc_env1, rc_env2_c)
		rcPredict1_c[i,] = predict.stm_point(spPosterior[i,cCols], rc_env1, rc_env2_c)
		rcPredict2_e[i,] = predict.stm_point(spPosterior[i,eCols], rc_env1_c, rc_env2)
		rcPredict2_c[i,] = predict.stm_point(spPosterior[i,cCols], rc_env1_c, rc_env2)

		# grid predictions
		grPredict_e[i,] = predict.stm_point(spPosterior[i,eCols], gr_env1, gr_env2)
		grPredict_c[i,] = predict.stm_point(spPosterior[i,cCols], gr_env1, gr_env2)
		grLambda[i,] = grPredict_c[i,] - grPredict_e[i,]
		grPres[i,] = as.integer(grLambda[i,] > 0)
	}
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
	saveRDS(respCurve, file.path('res','resp_curve',paste0(spName,'_respCurve.rds')))
	
	# maps
	spGrid = readRDS(file.path('res','sdm',paste0(spName, '_sdm_projection.rds')))
	spGrid$stm = apply(grLambda, 2, function(x) sum(x > 0, na.rm=TRUE)/length(x))
	spGrid$stm.var = spGrid$stm*(1-spGrid$stm) # binomial variance
	spGrid$sdmPres = as.integer(spGrid$sdm >= sdmThreshold)

	# categories for range disequilibrium
	spGrid$rde.present = apply(grPres, 2, function(x) sum(x == 1 & pGrid$sdmPres == 1)
	spGrid$rde.absent = apply(grPres, 2, function(x) sum(x == 0 & pGrid$sdmPres == 0)
	spGrid$rde.expand = apply(grPres, 2, function(x) sum(x == 1 & pGrid$sdmPres == 0)
	spGrid$rde.contract = apply(grPres, 2, function(x) sum(x == 0 & pGrid$sdmPres == 1)
	spGrid$rde = NA
	spGrid$rde[spGrid$rde.present >= spGrid$rde.expand & spGrid$rde.present >= spGrid$rde.contract] = 0
	spGrid$rde[spGrid$rde.expand > spGrid$rde.present & spGrid$rde.expand >= spGrid$rde.contract] = 1
	spGrid$rde[spGrid$rde.contract > spGrid$rde.present & spGrid$rde.contract >= spGrid$rde.expand] = 2
	saveRDS(spGrid, file.path('res','maps',paste0(spName,'_maps.rds')))
}