#!/usr/bin/Rscript

library(argparse)
library(rgdal)
library(fields)
library(reshape2)
library(raster)

# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
argList = parser$parse_args()
spName = argList$species

# set up projections
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers = CRS("+init=epsg:6350") # albers equal area conic NAD-83 north america
albers.lim = list(x=c(0, 3e6), y=c(0,3.5e6))



# get: rangeRaster, posteriorSummary
infile = paste("results/", spName, "/", spName, "_rangeRaster.rds", sep="")
rangeRaster = readRDS(infile)
proj4string(rangeRaster) = P4S.latlon
rangeRaster = projectRaster(rangeRaster, crs=P4S.albers)
# the transformation results in some impossible values, so we fix them
rangeRaster[rangeRaster < 0] = 0

infile = paste("results/", spName, "/", spName, "_posteriorSummary.rds", sep="")
posteriorSummary = readRDS(infile)
coordinates(posteriorSummary) = c('lon', 'lat')
gridded(posteriorSummary) = TRUE
posteriorRaster = brick(posteriorSummary)
setMinMax(posteriorRaster)
proj4string(posteriorRaster) = P4S.latlon
posteriorRaster = projectRaster(posteriorRaster, crs=P4S.albers)




# get some map data
ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, P4S.albers)
littleMap = readOGR(dsn=paste("dat/", spName, "/", spName, sep=""), layer=spName)
proj4string(littleMap) = P4S.latlon
littleMap = spTransform(littleMap, P4S.albers)
lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
# grab specific lakes
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
grLakes = spTransform(grLakes, P4S.albers)

plotbg = function(txt="", rangeMap=T)
{
	plot(ocean, col="white", add=T)
	plot(grLakes, col="white", add=T)
	mtext(txt, adj=0, cex = 1)
	
	if(rangeMap)
	plot(littleMap, border="red", add=T)

	# plot grid lines and labels
	llgridlines(ocean, easts=seq(-90, -40, 10), norths=seq(20,50,10), ndiscr=100, side="ES", offset=1e6)
	par(xpd=TRUE)
	westText = c("30°N", "40°N", "50°N")
	westCoords = list(rep(-3e5,3), c(0.75e6, 1.9e6, 3e6))
	southText = c("90°W", "80°W", "70°W")
	southCoords = list(c(0.65e6, 1.7e6, 2.8e6), rep(-0.14e6,3))
	cex.ll = 0.8
	text(westCoords[[1]], westCoords[[2]], westText, cex=cex.ll)
	text(southCoords[[1]], southCoords[[2]], southText, cex=cex.ll)
	par(xpd=FALSE)
}

rangeColors = colorRampPalette(c("#b30000", "#e34a33", "#fc8d59", "#fdcc8a", "#ffffff"), interpolate = 'spline', space='rgb', bias=1.0)(200)


w = 9
h = 3
pdf(width=w, height = h, file=paste("img/", spName, "/", spName, "_posterior_maps.pdf", sep=""))
## quartz(width=w, height = h)
par(mfrow=c(1,3), mar=c(2,2,2,1), oma=c(0,0,0,1))

# plot the posterior range boundary
## plot(rangeRaster, col=rev(rangeColors), main="pr(c - e == 0)", xaxt='n', yaxt='n')
## plotbg()

# plot the mean of c - e
lamColors = colorRampPalette(c("#008837", "#a6dba0", "#f7f7f7", "#c2a5cf", "#7b3294"), interpolate='spline', space="rgb", bias=1.0)(200)
# zrange = c(minValue(posteriorRaster$lambda), maxValue(posteriorRaster$lambda))
# zlims = round(zrange[which(abs(zrange) == min(abs(zrange)))] * c(-1, 1), 2)
# hist(posteriorRaster$lambda)

# zlims were chosen manually to make for pretty pictures
lam_zlims = list(
	"18032-ABI-BAL" = c(-0.15, 0.15),
	"19290-QUE-ALB" = c(-0.1, 0.1),
	"19481-BET-ALL" = c(-0.65, 0.65),
	"19489-BET-PAP" = c(-0.25, 0.25),
	"28728-ACE-RUB" = c(-0.2, 0.2),
	"28731-ACE-SAC" = c(-0.05, 0.05),
	"183302-PIC-MAR" = c(-0.15, 0.15),
	"183319-PIN-BAN" = c(-0.005, 0.005))
	
zlims = lam_zlims[[spName]]

plot(posteriorRaster$lambda, col=lamColors, zlim=zlims, main=paste(spName, ": c-e", sep=""), xaxt='n', yaxt='n')
plotbg()

lamMask = posteriorRaster$lambda < min(zlims)
presColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), interpolate='spline', bias=2, space="rgb")(200)
sdmColors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), interpolate='spline', bias=1, space="rgb")(200)

# plot the posterior uncertainty in c-e
plot(mask(posteriorRaster$lam.pres, lamMask, maskvalue=1), col=presColors, main="Proportion of samples where c > e", xaxt='n', yaxt='n')
plotbg()

# plot the sdm
plot(mask(posteriorRaster$sdm, lamMask, maskvalue=1), col=sdmColors, main="SDM", xaxt='n', yaxt='n')
plotbg()
dev.off()


