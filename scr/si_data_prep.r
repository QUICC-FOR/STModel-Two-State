#!/usr/bin/env Rscript

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/clim/plotClimate_unscaled.rds")
climScaling = readRDS("dat/clim/climate_scaling.rds")

# SDM data
for(spName in speciesList)
{
	spPath = file.path('dat', 'presence', paste0(spName, '_presence.rds'))
	spData = readRDS(spPath)
	if(!exists("rfData")) {
		rfData = spData
	} else
		rfData = merge(rfData, spData)
}
rfData = merge(rfData, climDat)[,-c(1,2)]
write.csv(rfData, "dat/si_rf_data.csv", row.names=FALSE)


# MCMC data
mcmcData = list()
for(spName in speciesList)
{
	mcmcDataFile = file.path('dat', 'mcmc', paste0(spName, '_mcmcDat.txt'))
	mcmcData[[spName]] = read.csv(mcmcDataFile)
	mcmcData[[spName]]$species = spName
	mcmcData[[spName]]$env1 = 
		(mcmcData[[spName]]$env1 * climScaling$scale['annual_mean_temp']) + 
		climScaling$center['annual_mean_temp']
	mcmcData[[spName]]$env2 = 
		(mcmcData[[spName]]$env2 * climScaling$scale['tot_annual_pp']) + 
		climScaling$center['tot_annual_pp']
}
mcmcData = do.call(rbind, mcmcData)
names(mcmcData)[6] = 'prevalence'
write.csv(mcmcData, "dat/si_mcmc_data.csv", row.names=FALSE)


# mortality & recruitment
demog = readRDS(file.path('res', 'demography', 'demography.rds'))
demog = demog[demog$interval <= 15,-c(1, 4)]
demog = demog[order(demog$species),]
write.csv(mcmcData, "dat/si_demog_data.csv", row.names=FALSE)

# species codes
speciesInfo = read.csv('dat/speciesInfo.csv')
speciesInfo = speciesInfo[which(speciesInfo$spName %in% speciesList), c('spName', 'genus', 'species')]
write.csv(speciesInfo, "dat/si_species.csv", row.names=FALSE)