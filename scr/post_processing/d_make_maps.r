#!/usr/bin/env Rscript

## depends:
##    res/posteriorGrid/*

## makes:
##    res/maps/*

library(coda)

speciesList = readRDS('dat/speciesList.rds')
models = c('0', 'g')

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
		
		spGrid$stm = spGrid$rde.present = spGrid$rde.expand = spGrid$rde.contract = NA
		for(i in 1:nrow(spGrid))
		{
			spGrid$stm[i] = sum(grLambda[,i] > 0, na.rm=TRUE)/length(grLambda[,i])
			spGrid$rde.present[i] = sum(grPres[,i] == 1 & spGrid$sdm.pres[i] == 1, na.rm=TRUE)/length(grPres[,i])
			spGrid$rde.expand[i] = sum(grPres[,i] == 1 & spGrid$sdm.pres[i] == 0, na.rm=TRUE)/length(grPres[,i])
			spGrid$rde.contract[i] = sum(grPres[,i] == 0 & spGrid$sdm.pres[i] == 1, na.rm=TRUE)/length(grPres[,i])
		}
		
		spGrid$stm.pres = as.integer(spGrid$stm > stmThreshold)
		spGrid$rde = (1 * (spGrid$stm.pres & spGrid$sdm.pres)) + 
				(2 * (spGrid$stm.pres & !spGrid$sdm.pres)) + 
				(3 * (!spGrid$stm.pres & spGrid$sdm.pres))
		# change absences to NA and subtract 1 from all cats
		spGrid$rde[spGrid$rde == 0] = NA 
		spGrid$rde = spGrid$rde - 1
		saveRDS(spGrid, file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
	}
}

