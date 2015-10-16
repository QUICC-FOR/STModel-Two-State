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

suppressWarnings(dir.create(file.path('img', 'mcmc'), recursive=TRUE))
suppressWarnings(dir.create(file.path('res', 'posterior'), recursive=TRUE))

remove_constant_cols = function(x, tol=1e-10)
{
	vars = apply(x, 2, var)
	cols = which(vars > tol)
	x[,cols]
}


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


# response curves
## do all species on one plot, draw a box around each pair of plots
## use this to accomplish the boxing:

## layout(matrix(c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16), 4, 4, byrow=TRUE))
## replicate(16, hist(rnorm(100)))
## par(xpd=NA)
## rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
##      grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
## rect( grconvertX(0.005, from='ndc'), grconvertY(0.005, from='ndc'),
##      grconvertX(0.495, from='ndc'), grconvertY(0.495, from='ndc'))
## rect( grconvertX(0.505, from='ndc'), grconvertY(0.505, from='ndc'),
##      grconvertX(0.995, from='ndc'), grconvertY(0.995, from='ndc'))
## rect( grconvertX(0.505, from='ndc'), grconvertY(0.005, from='ndc'),
##      grconvertX(0.995, from='ndc'), grconvertY(0.495, from='ndc'))


# STM posterior maps

# STM - SDM maps