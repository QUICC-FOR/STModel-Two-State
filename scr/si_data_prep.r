#!/usr/bin/env Rscript

# prepare data to be shared in the ms
library(randomForest)

dbInfo <- read.csv2('dat/db_info.csv', sep=';')
allowedDBs <- c('domtar_pp', 'us_pp')
allowedPlots <- dbInfo$plot_id[dbInfo$org_db_loc %in% allowedDBs]

speciesList = readRDS('dat/speciesList.rds')
climDat = readRDS("dat/clim/plotClimate_unscaled.rds")
climScaling = readRDS("dat/clim/climate_scaling.rds")

# SDM data; merged into one data frame, with climate
# also subset by plots we can publicly share
rfFiles <- sapply(speciesList, function(spName) file.path('dat', 'presence', paste0(spName, '_presence.rds')))
rfData <- Reduce(merge, lapply(rfFiles, readRDS))
rfData <- merge(rfData, climDat, all=FALSE)
rfData_public <- rfData[rfData$plot_id %in% allowedPlots,]
write.csv(rfData[,-c(1,2)], "dat/si_rf_data_ALL.csv", row.names=FALSE)
write.csv(rfData_public[,-c(1,2)], "dat/si_rf_data.csv", row.names=FALSE)


# MCMC data
mcmcFiles <- sapply(speciesList, function(spName) file.path('dat', 'stm_calib', paste0(spName, '_stm_calib.rds')))
stmCalib <- do.call(rbind, Map(function(path, name) cbind(readRDS(path), species=name), mcmcFiles, speciesList))
stmCalib$env1 <- (stmCalib$annual_mean_temp * climScaling$scale['annual_mean_temp']) + climScaling$center['annual_mean_temp']
stmCalib$env2 <- (stmCalib$tot_annual_pp * climScaling$scale['tot_annual_pp']) + climScaling$center['tot_annual_pp']
stmCalib$interval <- stmCalib$year2 - stmCalib$year1
stmCalib <- stmCalib[,c('plot', 'species', 'state1', 'state2', 'env1', 'env2', 'interval', 'prevalence')]
colnames(stmCalib) <- c('plot', 'species', 'initial', 'final', 'env1', 'env2', 'interval', 'prevalence')
stmCalib_public <- stmCalib[stmCalib$plot %in% allowedPlots, ]
write.csv(stmCalib[,-1], "dat/si_mcmc_data_ALL.csv", row.names=FALSE)
write.csv(stmCalib_public[,-1], "dat/si_mcmc_data.csv", row.names=FALSE)


# mortality & recruitment
demog = readRDS(file.path('res', 'demography', 'demography.rds'))
demog = demog[demog$interval <= 15,]
demog = demog[order(demog$species),]
demog_public <- demog[demog$plot_id %in% allowedPlots,]
write.csv(demog[,-c(1, 4)], "dat/si_demog_data_ALL.csv", row.names=FALSE)
write.csv(demog_public[,-c(1, 4)], "dat/si_demog_data.csv", row.names=FALSE)

# species codes
speciesInfo = read.csv('dat/speciesInfo.csv')
speciesInfo = speciesInfo[which(speciesInfo$spName %in% speciesList), c('spName', 'genus', 'species')]
write.csv(speciesInfo, "dat/si_species.csv", row.names=FALSE)
