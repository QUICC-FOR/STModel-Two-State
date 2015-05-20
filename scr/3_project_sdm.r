#!/usr/bin/Rscript

library(argparse)
library(gam)
library(randomForest)
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
argList = parser$parse_args()
spName = argList$species

load(paste("dat/", spName, "/", spName, "_processed.rdata", sep=""))
load(paste("results/", spName, "/", spName, "_sdm_models.rdata", sep=""))
climGrid = readRDS(paste("dat/", spName, "/", spName, "_climGrid_scaled.rds", sep=""))
# load(paste("dat/", spName, "/", spName, "_sdmData.rdata", sep=""))


# project the SDM prob of presence to the plots
transitionData.projected = transitionData.scaled
transitionData.projected$expectedGLM = predict(glm.mod, newdata = transitionData.scaled, type='response')
transitionData.projected$expectedGAM = predict(gam.mod, newdata = transitionData.scaled, type='response')
transitionData.projected$expectedRF = predict(rf.mod, newdata = transitionData.scaled, type='prob')[,2]
saveRDS(transitionData.projected, paste("dat/", spName, "/", spName, "_transitions_projected.rds", sep=""))


# project the SDMs to the climate grid
climGrid$glm.predict = predict(glm.mod, newdata=climGrid[,selectedVars], type='response')
climGrid$gam.predict = predict(gam.mod, newdata=climGrid[,selectedVars], type='response')
climGrid$rf.predict = predict(rf.mod, newdata=climGrid[,selectedVars], type='prob')[,2]
saveRDS(climGrid, paste("dat/", spName, "/", spName, "_climGrid_projected.rds", sep=""))






# make some plots

plotSettings = list(
	responseCurves = list(
		width = 12,
		height = 8,
		lwd = 1.5,
		col=list(glm='red', gam='blue', rf='cyan'),
		bty='n',
		mfrow = c(2, ceiling(length(selectedVars)/2))
	),
	map = list(
		zlim = c(0,1),
		sdmColors = colorRampPalette(c("#ffffff", "#ffffb2", "#fecc5c", "#fd8d3c", 
				"#e31a1c"), interpolate='spline', space="rgb", bias=1.3)(200),
		rangeBorder = '#FF3333bb',
		rangeLwd = 1.2,
		rangeCol = "#66666600",
		width = 12,
		height = 5.5
	)
)

## Response curves
with(plotSettings$responseCurves,
{
	if(interactive()) {
		quartz(w=width, h=height)
	} else {
		pdf(paste("img/", spName, "/", spName, "_response_curves.pdf", sep=""), w=width, h=height)
	}
	par(mfrow=mfrow)
	for(v in selectedVars) {
		newDat = as.data.frame(matrix(0, nrow=1000, ncol=length(selectedVars)))
		colnames(newDat) = selectedVars
		newDat[,v] = seq(-3,3, length.out=1000)

		plot(newDat[,v], predict(glm.mod, newdata=newDat, type='response'), lwd=lwd,
				type='l', col=col$glm, ylim=c(0,1), xlab=v, ylab="Pr(presence)", bty=bty)
		lines(newDat[,v], predict(gam.mod, newdata=newDat, type='response'), col=col$gam, 
				lwd=lwd)
		lines(newDat[,v], predict(rf.mod, newdata=newDat, type='prob')[,2], col=col$rf, 
				lwd=lwd)
		legend(0, 1, legend = c("GLM", "GAM", "RF"), lwd=1, col=c(col$glm, col$gam, 
				col$rf), bty='n')
	}
	if(!interactive()) dev.off()
})


# --------------------#
#  Make some maps
# --------------------#

library(rgdal)
library(fields)
with(plotSettings$map,
{
	ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
	lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
	# grab specific lakes
	lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
	grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]

	if(interactive()) {
		quartz(w=width, h=height)
	} else {
		pdf(paste("img/", spName, "/", spName, "_sdm_maps.pdf", sep=""), w=width, h=height)
	}

	layout(matrix(c(1,2,3,4,4,4), nrow=2, byrow=T), heights=c(1, 0.2))

	quilt.plot(climGrid$lon, climGrid$lat, climGrid$glm.predict, col=sdmColors, zlim=zlim, add.legend=F, xaxt='n', yaxt='n', useRaster=T, main="GLM")
		plot(ocean, col="white", add=T)
		plot(grLakes, col="white", add=T)
	quilt.plot(climGrid$lon, climGrid$lat, climGrid$gam.predict, col=sdmColors, zlim=zlim, add.legend=F, xaxt='n', yaxt='n', useRaster=T, main="GAM")
		plot(ocean, col="white", add=T)
		plot(grLakes, col="white", add=T)
	quilt.plot(climGrid$lon, climGrid$lat, climGrid$rf.predict, col=sdmColors, zlim=zlim, add.legend=F, xaxt='n', yaxt='n', useRaster=T, main="RF")
		plot(ocean, col="white", add=T)
		plot(grLakes, col="white", add=T)
	par(mar=c(3,3,3,3))
	image(x=seq(zlim[1],zlim[2],length.out=101), z=matrix(seq(zlim[1], zlim[2],length.out=100), 
		nrow=100, ncol=1), zlim=zlim, col=sdmColors, yaxt='n', xlab='', ylab='', useRaster=T)
	mtext("Probability of Presence", line=0.5, cex = 0.7)
	if(!interactive()) dev.off()
})

