#!/usr/bin/env Rscript
# compute the posterior distribution of the response curves

## depends:
##    res/posterior/*

## makes:
##    res/resp_surface.rds

library(coda)
library(parallel)

# number of points in each dimension of response surface
# number of posterior rows to use
rsRes = 100
posterior.n = 1000

speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
source('scr/stm_functions.r')

climScale = readRDS('dat/clim/climate_scaling.rds')
climDat = readRDS('dat/clim/plotClimate_scaled.rds')

# the posterior columns containing extinction and colonization parameters
e_pars = function(pars)	if(length(pars) > 2) pars[6:10] else pars[2]
c_pars = function(pars)	if(length(pars) > 2) pars[1:5] else pars[1]

# set scales for the response surface
rs.env1.o = seq(-5, 25, length.out=rsRes) # temperature range on original scale
rs.env2.o <- seq(500, 2000, length.out=rsRes) # precip range on original scale
rs.env1 <- scale(rs.env1.o, center = climScale$center['annual_mean_temp'], scale = climScale$scale['annual_mean_temp'])
rs.env2 = scale(rs.env2.o, center = climScale$center['tot_annual_pp'], scale = climScale$scale['tot_annual_pp'])
rsDat <- expand.grid(annual_mean_temp=rs.env1, tot_annual_pp=rs.env2)

respSurfaceList <- mclapply(speciesList, function(spName) {
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_0_samples.rds'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]
	info = speciesInfo[speciesInfo$spName == spName,]
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'_stm_calib.rds')))

	rsPredict <- rowMeans(apply(samples, 1, function(pars) {
		predict.stm_point(c_pars(pars), rsDat[,1], rsDat[,2]) - predict.stm_point(e_pars(pars), rsDat[,1], rsDat[,2])
	}))

	matrix(rsPredict, nrow=length(rs.env1), ncol=length(rs.env2))
})
names(respSurfaceList) <- speciesList

saveRDS(list(x.precip = rs.env2.o, y.temp = rs.env1.o, z =  respSurfaceList), "res/resp_surface.rds")

