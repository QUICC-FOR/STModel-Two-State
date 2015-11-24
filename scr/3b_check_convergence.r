#!/usr/bin/env Rscript

## run this to check the convergence of chains; prints the gelman diagnostic to the screen

library(coda)
library(foreach)

# full model species -- edit to use whatever species are desired
spList.full = c('505490-THU-OCC', '183319-PIN-BAN', '183412-LAR-LAR')

# int only species
spList.int = c('19287-QUE-MAC', '183397-TSU-CAN')


get.samples = function(sp, mod)
{
	chains = c('ch1', 'ch2', 'ch3')
	foreach(ch = chains, .final=mcmc.list) %do%
	{
		modname = paste0(mod, '_', ch)	
		samples = read.csv(file.path('res','mcmc',sp,modname,'posterior.csv'))
		constcols = which(apply(samples, 2, function(x) length(unique(x)) == 1))
		mcmc(samples[,-constcols])
	}
}

samples.full = foreach(sp = spList.full, .final=function(x) { names(x) = spList.full; x }) %do% {
	get.samples(sp, '0')
}

lapply(samples.full, gelman.diag)


samples.int = foreach(sp = spList.int, .final=function(x) { names(x) = spList.int; x }) %do% {
	get.samples(sp, '0i')
}

lapply(samples.int, gelman.diag)