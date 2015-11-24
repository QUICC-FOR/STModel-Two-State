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

suppressWarnings(dir.create(file.path('res', 'posterior'), recursive=TRUE))


remove_constant_cols = function(x, tol=1e-10)
{
	vars = apply(x, 2, var)
	cols = which(vars > tol)
	x[,cols]
}

plot_mcmc_summary = function(PD, base.name, make.png = TRUE)
{
	dens = FALSE
	tr = TRUE
	w = 10
	h = 15
	dpi=600
	fontsize=11
	fname = file.path('img','mcmc', paste0(base.name, '_traceplots.png'))
	
	if(make.png)
	{
		png(width=as.integer(dpi*w), height=as.integer(dpi*h), file=fname, 
				pointsize=fontsize, res=dpi)
	} else
	{
		quartz(width = w, height=h)
	}

	if(ncol(PD) == 2)
	{
		pl.dim = c(2,2)
	} else {
		pl.dim = c(5,4)
	}
	
	par(mfrow=pl.dim, mar=c(3,3,0,0), mgp=c(1.5, .5, 0), oma=c(0,0,2,2), cex=0.5, bty='n')
	plot(PD, ask=FALSE, auto.layout=FALSE, main='')
	if(make.png) dev.off()
}


for(spName in speciesList)
{
	cat(paste("Reading posterior for", spName, "\n"))
	info = speciesInfo[speciesInfo$spName == spName,]
	models = c('0', 'i0')
	thin = c('0' = info$thin, 'i0'=5)
	burnin = c('0' = info$burnin, 'i0'=10000)
	for(mod in models)
	{
		base.name = paste(spName, mod, sep='_')
		samples.name = file.path('res', 'posterior', paste0(base.name, '_samples.rds'))
		mod.path = file.path('res','mcmc',spName,mod,'posterior.csv')
		samples = read.csv(mod.path)
		samples = remove_constant_cols(samples)
		samp.mcmc = mcmc(samples, start=thin[[mod]]*burnin[[mod]]+1, thin=thin[[mod]])
		saveRDS(samp.mcmc, samples.name)
		plot_mcmc_summary(samp.mcmc, base.name)
	}
}
