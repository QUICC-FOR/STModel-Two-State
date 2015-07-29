#!/usr/bin/Rscript
library(coda)
library(raster)

## parameters = commandArgs(trailingOnly = TRUE)
## spName = parameters[1]
## if(length(parameters) > 1) {
## 	curvePrecision = parameters[2]
## } else {
## 	curvePrecision = 500
## }

spList = c('19231-CAR-GLA', '19277-QUE-FAL', '183397-TSU-CAN', '18048-JUN-VIR', '505490-THU-OCC')
speciesInfo = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE, colClasses='character')
load('dat/map_projections.rdata')

compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}

for(spName in spList)
{
	cat(paste("\nStarting", spName, "\n"))
	spInfo = speciesInfo[speciesInfo$spName == spName,]
	env1 = spInfo$env1
	env2 = spInfo$env2
	design = spInfo$design
	design = sapply(1:nchar(design), function(i) as.integer(substr(design, i, i)))



	posterior = readRDS(file.path('species', spName, 'res', paste(spName, 'posterior_thinned.rds', sep='_')))
	climGrid = readRDS('dat/climateGrid_scaled.rds')
	parBase = rep(0, length(design))

	x1 = seq(min(climGrid[,env1]), max(climGrid[,env1]), length.out=curvePrecision)
	x2 = seq(min(climGrid[,env2]), max(climGrid[,env2]), length.out=curvePrecision)

	# now compute the y values for each x value for every posterior sample
	# first we have to choose a value at which to fix the second environmental variable
	# do this with the mean values of the parameters
	meanParams = parBase
	meanParams[design == 1] = summary(posterior)$statistics[,1]

	# find the point at which env1 maximizes lambda
	env1.yc = compute_c(meanParams, x1, rep(0, length(x1)))
	env1.ye = compute_e(meanParams, x1, rep(0, length(x1)))
	env1.lam = env1.yc - env1.ye
	env1.maxx = x1[which(env1.lam == max(env1.lam))[1]]

	# find the point at which env2 maximizes lambda with the max from the prev step
	env2.yc = compute_c(meanParams, rep(env1.maxx, length(x2)), x2)
	env2.ye = compute_e(meanParams, rep(env1.maxx, length(x2)), x2)
	env2.lam = env2.yc - env2.ye
	env2.maxx = x2[which(env2.lam == max(env2.lam))[1]]

	# re-do the computation of x1 with the new max from x2
	env1.yc = compute_c(meanParams, x1, rep(env2.maxx, length(x1)))
	env1.ye = compute_e(meanParams, x1, rep(env2.maxx, length(x1)))
	env1.lam = env1.yc - env1.ye
	env1.maxx = x1[which(env1.lam == max(env1.lam))[1]]

	y1.c = t(do.call(cbind, lapply(posterior, function(x)
	{
		sapply(1:nrow(x), function(i)
		{
			params = parBase
			params[design == 1] = x[i,]
			compute_c(params, x1, rep(env2.maxx, length(x1)))
		})
	})))


	y1.e = t(do.call(cbind, lapply(posterior, function(x)
	{
		sapply(1:nrow(x), function(i)
		{
			params = parBase
			params[design == 1] = x[i,]
			compute_e(params, x1, rep(env2.maxx, length(x1)))
		})
	})))


	y2.c = t(do.call(cbind, lapply(posterior, function(x)
	{
		sapply(1:nrow(x), function(i)
		{
			params = parBase
			params[design == 1] = x[i,]
			compute_c(params, rep(env1.maxx, length(x2)), x2)
		})
	})))


	y2.e = t(do.call(cbind, lapply(posterior, function(x)
	{
		sapply(1:nrow(x), function(i)
		{
			params = parBase
			params[design == 1] = x[i,]
			compute_e(params, rep(env1.maxx, length(x2)), x2)
		})
	})))

	# now put the x variables back on their original scales
	climScale = readRDS("dat/climate_scaling.rds")
	x1.us = (x1 * climScale$scale[env1]) + climScale$center[env1]
	x2.us = (x2 * climScale$scale[env2]) + climScale$center[env2]

	responseCurves = data.frame(
		env1 = x1.us,
		env2 = x2.us,
		col1.mean = colMeans(y1.c),
		col1.lo = apply(y1.c, 2, quantile, 0.025),
		col1.up = apply(y1.c, 2, quantile, 0.975),
		ext1.mean = colMeans(y1.e),
		ext1.lo = apply(y1.e, 2, quantile, 0.025),
		ext1.up = apply(y1.e, 2, quantile, 0.975),
		col2.mean = colMeans(y2.c),
		col2.lo = apply(y2.c, 2, quantile, 0.025),
		col2.up = apply(y2.c, 2, quantile, 0.975),
		ext2.mean = colMeans(y2.e), 
		ext2.lo = apply(y2.e, 2, quantile, 0.025),
		ext2.up = apply(y2.e, 2, quantile, 0.975))
	
	saveRDS(responseCurves, file.path('species', spName, 'res', paste(spName, 'responseCurves.rds', sep='_')))




	## now do the map data


	prj.ras = function(x)
	{
		coordinates(x) = 1:2
		gridded(x) = TRUE
		x = raster(x)
		proj4string(x) = P4S.latlon
		projectRaster(x, crs=stmMapProjection)
	}


	# these are quite heavy, so they have to be done one grid cell at a time to conserve memory
	# the intermediate data structures are much larger than the result

	map_data = function(post.list, e1, e2)
	{
		# note that pos.list should be an mcmc list with equal-sized chunks of parallel chains
		# e1 and e2 should be scalars - we are doing this just for a single point on the grid
	
		pos = do.call(rbind, post.list)
		md = as.data.frame(t(sapply(1:nrow(pos), function(i)
		{
				params = parBase
				params[design == 1] = pos[i,]
				c(c=compute_c(params, e1, e2),
					e=compute_e(params, e1, e2))
		})))

		md$lam = md$c - md$e

		c(c=mean(md$c), e=mean(md$e), lam = mean(md$lam), pres = sum(md$lam > 0)/nrow(md))
	}

	mapData = matrix(NA, nrow=nrow(climGrid), ncol = 6)
	mapData[,1] = climGrid$lon
	mapData[,2] = climGrid$lat
	for(i in 1:nrow(climGrid))
	{
		mapData[i,3:6] = map_data(posterior, climGrid[i,env1], climGrid[i,env2])
		if(i %% 1000 == 0) cat(paste("  ", Sys.time(), "-- Finished", i, "of", nrow(climGrid), "map cells\n"))
	}


	mapData.c = data.frame(lon = mapData[,1], lat = mapData[,2], c = mapData[,3])
	mapData.e = data.frame(lon = mapData[,1], lat = mapData[,2], e = mapData[,4])
	mapData.lam = data.frame(lon = mapData[,1], lat = mapData[,2], lam = mapData[,5])
	mapData.pres = data.frame(lon = mapData[,1], lat = mapData[,2], pres = mapData[,6])
	mapData.sdm = readRDS(file.path('species', spName, 'res', paste(spName, 'sdm_grid_projection.rds', sep='_')))
		
	grid.c = prj.ras(mapData.c)
	grid.e = prj.ras(mapData.e)
	grid.lam = prj.ras(mapData.lam)
	grid.pres = prj.ras(mapData.pres)
	grid.sdm = prj.ras(mapData.sdm)

	save(grid.c, grid.e, grid.lam, grid.pres, grid.sdm, 
			file=file.path('species', spName, 'res', paste(spName, 'mapRasters.rdata', sep='_')))
}
