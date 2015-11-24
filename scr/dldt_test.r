#! /usr/bin/env Rscript

library(numDeriv)
library(rootSolve)
library(foreach)
library(iterators)
library(coda)
library(doParallel)
registerDoParallel(cores=detectCores())

## library(ellipse)
speciesInfo = read.csv('dat/speciesInfo.csv')
climScale = readRDS('dat/clim/climate_scaling.rds')
dldt.path = file.path('res','dldt.rds')
speciesList = readRDS('dat/speciesList.rds')
posterior.n = 1000
climGrid = readRDS(file.path('dat', 'clim', 'climateGrid_scaled.rds'))

recompute = FALSE
arg = commandArgs(trailingOnly = TRUE)
if('--recompute' %in% arg | '-r' %in% arg) recompute = TRUE


lam_t = function(T, pars, P)
{
	Cl = pars[1] + pars[2]*T + pars[3]*P + pars[4]*T^2 + pars[5]*P^2
	El = pars[6] + pars[7]*T + pars[8]*P + pars[9]*T^2 + pars[10]*P^2
	plogis(Cl) - plogis(El)
}


get_derivs = function(curPars, temp.range, precip, stepSize = 0.25)
{
	roots = uniroot.all(lam_t, pars=curPars, P = precip, lower=temp.range[1], upper=temp.range[2])
	if(length(roots) > 0)
	{
		dldt = grad(lam_t, roots, pars = curPars, P = precip)
		type = sapply(roots, function(r)
			# if we get colder (root - 0.01) and lambda is negative, then it is a northern limit
			ifelse(lam_t(r-0.01, curPars, precip) < 0, "northern", "southern"))
	} else {
		dldt = roots = type = NA
	}

	data.frame(dldt=dldt, root=roots, type=type, stringsAsFactors = FALSE)
}


compute_species_dldt = function(spName)
{
	# get the median precipitation experienced at the range boundary
	# we use the 'uncertainty range' where p is between 0.05 and 0.95
	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	spGrid = merge(spGrid, climGrid[,c('x','y', 'annual_mean_temp', 'tot_annual_pp')])
	precip = spGrid$tot_annual_pp[spGrid$stm>0.05 & spGrid$stm < 0.95]
	pp.range = seq(min(precip, na.rm=TRUE), max(precip, na.rm=TRUE), length.out=100)
	trange = range(spGrid$annual_mean_temp)
	
	# get posterior
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_0_samples.rds'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]

	info = speciesInfo[speciesInfo$spName == spName,]
	cat("starting pp iteration from", pp.range[1], "to", pp.range[length(pp.range)], "\n")
	foreach(pp = pp.range, .combine=rbind) %dopar%
	{	
		vals = foreach(curPars=iter(samples, by='row'), .combine=rbind) %do%
				get_derivs(curPars, trange, pp)
		vals = vals[complete.cases(vals),]
		vals$precip = pp
		vals
	}	
}

spName = "505490-THU-OCC"
dldt.df = compute_species_dldt(spName)
saveRDS(dldt.df, "res/dldt_test.rds")

dldt.df$dldt = abs(dldt.df$dldt)

dldt.mean = tapply(dldt.df$dldt, dldt.df$type, mean)
temp.mean = tapply(dldt.df$root, dldt.df$type, mean)
dldt.lower = tapply(dldt.df$dldt, dldt.df$type, quantile, 0.05)
dldt.upper = tapply(dldt.df$dldt, dldt.df$type, quantile, 0.95)
temp.lower = tapply(dldt.df$root, dldt.df$type, quantile, 0.05)
temp.upper = tapply(dldt.df$root, dldt.df$type, quantile, 0.95)

pdf(file="dldt_test.pdf")
plot(temp.mean, abs(dldt.mean), xlim=c(-2, 5), ylim=c(0,0.8), col=c('blue', 'red'), pch=16)
segments(temp.mean, dldt.lower, temp.mean, dldt.upper)
segments(temp.lower, abs(dldt.mean), temp.upper, abs(dldt.mean))
dev.off()