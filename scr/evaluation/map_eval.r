#!/usr/bin/env Rscript

# gets posterior tss/roc at the map scale for a single species/model

# depends:
#    res/posterior/*
#    dat/eval/*

# produces:
#    res/eval/SPNAME_MOD_posterior_eval.rds


library(foreach)
library(doParallel)
library(biomod2)

clArgs <- commandArgs(trailingOnly = TRUE)

## these will be read in from the command line
spName = clArgs[1]
modName = if(length(clArgs) > 1) clArgs[2] else '0'
numCores = if(length(clArgs) > 2) clArgs[3] else detectCores()
sampleSize = if(length(clArgs) > 3) clArgs[4] else NA

source('scr/stm_functions.r')
registerDoParallel(cores=numCores)

cat("Evaluation for ", spName, " model ", modName, "\n\n")
cat("Evaluation for ", spName, " model ", modName, "\n\n", file = stderr())

# get STM fit
posterior = readRDS(file.path('res', 'posterior', paste0(spName, '_posterior.rds')))[[modName]]
if(is.null(sampleSize) || is.na(sampleSize)) {
	samples = posterior
} else
	samples = posterior[sample(nrow(posterior), sampleSize),]


mapValidCells = readRDS(file.path('dat', 'eval', paste0(spName, '_evalCells.rds')))
env1 = mapValidCells$annual_mean_temp
env2 = mapValidCells$tot_annual_pp

eval.posterior <- foreach(pars = iter(samples, by='row'), .combine=rbind, 
		.final = function(x) {
			data.frame(tss = x[,1], roc = x[,2])
}) %dopar% {
	if(length(pars) == 2)
	{
		cp = pars[1]
		ep = pars[2]
	} else {
		cp = pars[1:5]
		ep = pars[6:10]
	}
	C = predict.stm_point(cp, env1, env2)
	E = predict.stm_point(ep, env1, env2)
	fit = as.integer(C > E)
	c(Find.Optim.Stat(Stat="TSS", Fit=(1000*fit), Obs=mapValidCells$obs)[1],
	Find.Optim.Stat(Stat="ROC", Fit=(1000*fit), Obs=mapValidCells$obs)[1])
}
saveRDS(eval.posterior, file.path('res', 'eval', paste0(spName, '_', modName, '_posterior_eval.rds')))