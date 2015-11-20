#!/usr/bin/env Rscript

## reads all posterior results in and processes them into CODA objects; also makes some plots
## depends:
##    res/posteriorGrid/*

## makes:
##    res/maps/*

library(coda)
library(foreach)
library(doParallel)
numCores=detectCores() - 2

speciesList = readRDS('dat/speciesList.rds')
models = c('0', 'i0', 'g', 'ig')

arg = commandArgs(TRUE)
if(length(arg) > 0)
{
	sps = which(speciesList %in% arg)
	if(length(sps) == 0) 
	{
		warning("No species specified on command line; falling back to default")
	} else 
		speciesList = speciesList[which(speciesList %in% arg)]
	corTest = as.numeric(arg[length(arg)])
	if(!is.na(corTest)) numCores = corTest
}
registerDoParallel(cores=numCores)

if(length(speciesList) == 0) stop("Error: no species specified")
cat("Will process posteriors for:\n")
for(sp in speciesList) cat("  ", sp, "\n")

suppressWarnings(dir.create(file.path('res', 'maps'), recursive=TRUE))
sdmThreshold = 0.1
stmThreshold = 0.1




for(spName in speciesList)
{
	for(mod in models)
	{
		postGrid = readRDS(file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
		spGrid = readRDS(file.path('res','sdm',paste0(spName, '_sdm_projection.rds')))
		spGrid$sdm.pres = as.integer(spGrid$sdm >= sdmThreshold)
		grPres = postGrid$stmPres
		grLambda = postGrid$lambda
		
		spGrResult = foreach(i = 1:nrow(spGrid), .combine=rbind, .final=data.frame) %dopar% {
			c(stm = sum(grLambda[,i] > 0, na.rm=TRUE)/length(grLambda[,i]),
				rde.present = sum(grPres[,i] == 1 & spGrid$sdm.pres[i] == 1, na.rm=TRUE)/length(grPres[,i]),
				rde.expand = sum(grPres[,i] == 1 & spGrid$sdm.pres[i] == 0, na.rm=TRUE)/length(grPres[,i]),
				rde.contract = sum(grPres[,i] == 0 & spGrid$sdm.pres[i] == 1, na.rm=TRUE)/length(grPres[,i])
			)}
		spGrResult$stm.pres = as.integer(spGrResult$stm > stmThreshold)
		spGrResult$rde = (1 * (spGrResult$stm.pres & spGrid$sdm.pres)) + 
				(2 * (spGrResult$stm.pres & !spGrid$sdm.pres)) + 
				(3 * (!spGrResult$stm.pres & spGrid$sdm.pres))
		# change absences to NA and subtract 1 from all cats
		spGrResult$rde[spGrResult$rde == 0] = NA 
		spGrResult$rde = spGrResult$rde - 1
		saveRDS(spGrResult, file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
	}
}

