#!/usr/bin/RScript

library(argparse)
library(coda)
library(sp)
library(raster)
library(reshape2)

# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-r", "--raster", default=0, type="integer", help="raster resolution (default 0)")
argList = parser$parse_args()
spName = argList$species

outDir = paste("results/", spName, "/", sep="")
rasterDim = as.integer(argList$raster)


timestamp = function() paste(Sys.time(), spName, "--")

compute_e = function(p, env1, env2)
{
	plogis(p[6] + env1*p[7] + env2*p[8] + env1^2*p[9] + env2^2*p[10])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5])
}

compute_boundary = function(dat, curCount, slow=FALSE)
{
	dat_wide = acast(dat,lon ~ lat, value.var = "lambda")
	xc = contourLines(x = sort(unique(dat$lon)), y = sort(unique(dat$lat)), z = dat_wide, levels=0)
	
	# sometimes there is no range boundary (either it is extinct everywhere or occurs everywhere)
	if(length(xc) == 0)
	return(curCount)

	if(!slow)
	{
		xcul = lapply(xc, function(x) as.matrix(cbind(x$x, x$y)))
		xcuul = as.data.frame(do.call(rbind, xcul))
		coordinates(xcuul) = c(1,2)
		boundary = rasterize(xcuul, curCount)
	} else {
		# this is more correct, but takes 10x as long
		xcul = lapply(xc, function(x) data.frame(x=x$x, y=x$y))
		xLine = lapply(xcul, Line)
		xLines = list(Lines(xLine, "1"))
		xSLines = SpatialLines(xLines)
		boundary = rasterize(xSLines, curCount)
	}
	vl = values(boundary)
	vl[is.na(vl)] = 0
	vl[vl != 0] = 1
	values(boundary) = vl
	return(boundary + curCount)
}



cat(paste("\n", timestamp(), " Starting posterior processing for species ", spName, "\n", sep = ""))

# read data
cat(paste(timestamp(), "converting raw MCMC data into CODA object\n"))
postFile = paste(outDir, "posterior.csv", sep="")
posterior = read.table(postFile, sep=',', header=TRUE, row.names=NULL)
posterior = mcmc(posterior, thin=25, start = 25*5000)
saveRDS(posterior, file=paste(outDir, spName, "_posterior.rds", sep=""))

climName = paste("dat/", spName, "/", spName, "_climGrid_projected.rds", sep="")
climDat = readRDS(climName)
if(spName == "183319-PIN-BAN")
{
	sdm = climDat$gam.predict
} else 
{
	sdm = climDat$rf.predict
}

if(rasterDim > 0)
{
	# set up a raster to hold the range boundary
	cat(paste(timestamp(), "Setting up raster\n"))
	marginGrid = data.frame(
		lon = climDat$lon,
		lat = climDat$lat,
		count = 0)
	coordinates(marginGrid) = ~lon+lat
	gridded(marginGrid) = TRUE
	rangeRaster = raster(marginGrid)
	if(rasterDim > 1)
	{
		rangeRaster = aggregate(rangeRaster, fact=as.integer(rasterDim), fun = sum)
	}
}

cat(paste(timestamp(), "computing posterior statistics on", nrow(posterior), "samples\n"))
posteriorData = list()
posteriorData$lambda = posteriorData$other = matrix(nrow = nrow(climDat), ncol=nrow(posterior))
for(i in 1:nrow(posterior))
{
	er = compute_e(posterior[i,], climDat$annual_mean_temp,  climDat$tot_annual_pp)
	cr = compute_c(posterior[i,], climDat$annual_mean_temp,  climDat$tot_annual_pp)
	posteriorData$lambda[,i] = cr - er
	posteriorData$other[,i] = cr*sdm*(1 - sdm) - er*sdm


	if(rasterDim > 0)
	{
		rangeRaster = compute_boundary(data.frame(lon = climDat$lon, lat = climDat$lat, 
			lambda = posteriorData$lambda[,i]), rangeRaster)
	}

	if(i %% 500 == 0)
	{
		cat(paste("    ", timestamp(), "completed", i, "of", nrow(posterior), "\n"))
		flush.console()
	}
}

posteriorData$pres = (posteriorData$lambda > 0)
if(rasterDim > 0) rangeRaster = rangeRaster / nrow(posterior)
saveRDS(rangeRaster, file = paste(outDir, spName, "_rangeRaster.rds", sep=""))



cat(paste(timestamp(), "computing posterior summary statistics\n"))
posteriorSummary = data.frame(
	lon = climDat$lon,
	lat = climDat$lat,
	e = compute_e(colMeans(posterior), climDat$annual_mean_temp, climDat$tot_annual_pp),
	c = compute_c(colMeans(posterior), climDat$annual_mean_temp, climDat$tot_annual_pp),
	sdm = sdm,
	lam.sd = apply(posteriorData$lambda, 1, sd),
	lam.cv = apply(posteriorData$lambda, 1, function(x) abs(sd(x)/mean(x))),
	lam.pres = rowSums(posteriorData$pres)/nrow(posterior),
	oth.sd = apply(posteriorData$other, 1, sd),
	oth.cv = apply(posteriorData$other, 1, function(x) abs(sd(x)/mean(x))))
posteriorSummary = within(posteriorSummary,	
{
	lambda <- c-e
	other <- c*sdm*(1-sdm) - e*sdm    ## what should this parameter be called??
})

outPS = paste(outDir, spName, "_posteriorSummary.rds", sep="")
saveRDS(posteriorSummary, file=outPS)

cat(paste(timestamp(), "done\n"))
