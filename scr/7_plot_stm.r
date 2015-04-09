library(argparse)
library(fields)
library(gam)

parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-c", "--cubic", action="store_true", default=FALSE, help="use cubic terms in model (default is square only)")
parser$add_argument("-r", "--rf", action="store_true", default=FALSE, help="random forest for prevalence (default: GLM)")
parser$add_argument("-g", "--gam", action="store_true", default=FALSE, help="use GAM for prevalence (default: GAM)")
parser$add_argument("-f", "--fraction", default=1, type="double", help="proportion of data to use for model fitting")
argList = parser$parse_args()
argList$species = "18032-ABI-BAL"


annealResults = readRDS(paste("results/", argList$species, "_anneal_parameters_", argList$fraction, ".rds", sep=""))
load(paste("dat/", argList$species, "_sdmData.rdata", sep=""))
load(paste("results/", argList$species, "_sdm_models.rdata", sep=""))
climGrid = readRDS(paste("dat/", argList$species, "_climGrid_projected.rds", sep=''))
#climGrid = readRDS("dat/climPast_scaled.rds")



plotSettings = list(
	responseCurve = list(
		width = 12,
		height = 6,
		lwd = 1.5,
		col=list(gamma='blue', epsilon = 'red', sdm = 'black', col = 'cyan'),
		bty='n'
	),
	map = list(
		zlim = c(0,1),
		rangeBorder = '#FF3333bb',
		rangeLwd = 1.2,
		rangeCol = "#66666600",
		width = 12,
		height = 6,
		ceCols = c('white', '#2FB59A')
	)
)




p = annealResults$par


with(plotSettings$responseCurve,
{
	if(interactive()) {
		quartz(w=width, h=height)
	} else {
		pdf(paste("img/", argList$species, "_stm_response_curves.pdf", sep=""), w=width, h=height)
	}
	par(mfrow=c(1,2))
	
	stm_newdat = function(v, variables)
	{
		newDat = as.data.frame(matrix(0, nrow=1000, ncol=length(variables)))
		colnames(newDat) = variables
		newDat[,v] = seq(-3,3, length.out=1000)
		return(newDat)
	}
	
	stm_project = function(v, pars, newDat)
	{

		probs = data.frame(
			col = plogis(pars[1] + pars[2]*newDat[,v] + pars[3] * newDat[,v]^2),
			ext = plogis(pars[4] + pars[5]*newDat[,v] + pars[6] * newDat[,v]^2),
			pres = predict(glm.mod, newdata=newDat, type='response'))

		if(argList$gam) {
			require(gam)
			probs$pres = predict(gam.mod, newdata=newDat, type='response')
		} else if(argList$rf) {
			require(randomForest)
			probs$pres = predict(rf.mod, newdata=newDat, type='prob')[,2]
		}

		if(argList$cubic)
		{
			probs = within(probs, 
			{
				col = plogis(p[1] + p[2]*newDat[,v] + p[3]*newDat[,v]^2 + p[4]*newDat[,v]^3)
				ext = plogis(p[5] + p[6]*newDat[,v] + p[7]*newDat[,v]^2 + p[8]*newDat[,v]^3)
			})
		}
		return(probs)
	}

	newDat = stm_newdat('annual_mean_temp', selectedVars)
	probs = stm_project('annual_mean_temp', p[c(1, 2, 4, 6, 7, 9)], newDat)
	plot(newDat$annual_mean_temp, probs$col, type='l', col=col$gamma, lwd=lwd, bty=bty, xlab="Temperature", ylab="Prob", ylim=c(1e-10,1), log="y")
	lines(newDat$annual_mean_temp, probs$ext, col=col$epsilon, lwd=lwd)
	lines(newDat$annual_mean_temp, probs$pres, col=col$sdm, lwd=lwd*0.5)
	lines(newDat$annual_mean_temp, probs$pres*probs$col, col=col$col, lwd=lwd)
	legend(1, 1e-7, c(expression(gamma), expression(epsilon), "SDM", expression(gamma %*% SDM)), lwd=c(1,1,.5, 1), col=c(col$gamma, col$epsilon, col$sdm, col$col))

	newDat = stm_newdat('tot_annual_pp', selectedVars)
	probs = stm_project('tot_annual_pp', p[c(1, 3, 5, 6, 8, 10)], newDat)
	plot(newDat$tot_annual_pp, probs$col, type='l', col=col$gamma, lwd=lwd, bty=bty, xlab="Precipitation", ylab="Prob", ylim=c(1e-10,1), log="y")
	lines(newDat$tot_annual_pp, probs$ext, col=col$epsilon, lwd=lwd)
	lines(newDat$tot_annual_pp, probs$pres, col=col$sdm, lwd=lwd*0.5)
	lines(newDat$tot_annual_pp, probs$pres*probs$col, col=col$col, lwd=lwd)
})


with(plotSettings$map,
{
	if(interactive()) {
		quartz(w=width, h=height)
	} else {
		pdf(paste("img/", argList$species, "_stm_presence_maps.pdf", sep=""), w=width, h=height)
	}
	par(mfrow=c(1,2))

	library(rgdal)
	ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
	lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
	mapleRange = readOGR(dsn="dat/acersacr", layer="acersacr")
	# grab specific lakes
	lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
	grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
	
	climGrid$sdm = predict(gam.mod, newdata=climGrid[,selectedVars], type='response')
	climGrid$gamma = plogis(p[1] + p[2]*climGrid$annual_mean_temp + 
			p[3]*climGrid$tot_annual_pp + p[4]*climGrid$annual_mean_temp^2 + 
			p[5]*climGrid$tot_annual_pp^2)
	climGrid$epsilon = plogis(p[6] + p[7]*climGrid$annual_mean_temp + 
			p[8]*climGrid$tot_annual_pp + p[9]*climGrid$annual_mean_temp^2 + 
			p[10]*climGrid$tot_annual_pp^2)
	climGrid$colonization = climGrid$sdm * climGrid$gamma
	climGrid$diff.gamma = climGrid$gamma - climGrid$epsilon
	climGrid$diff.col = climGrid$colonization - climGrid$epsilon
	climGrid$status = as.integer(ifelse(climGrid$diff.col <= 0, 0, 1))
	climGrid$status.gamma = as.integer(ifelse(climGrid$diff.gamma <= 0, 0, 1))

	quilt.plot(climGrid$lon, climGrid$lat, climGrid$status, col=ceCols, add.legend=F, xaxt='n', yaxt='n', useRaster=T, main="Prevalence = SDM")
		plot(ocean, col="white", add=T)
		plot(mapleRange[c(1,3,47,43),], border=rangeBorder, col=rangeCol, add=T, lwd=rangeLwd)
		plot(grLakes, col="white", add=T)
		
	quilt.plot(climGrid$lon, climGrid$lat, climGrid$status.gamma, col=ceCols, add.legend=F, xaxt='n', yaxt='n', useRaster=T, main="Prevalence = 1")
		plot(ocean, col="white", add=T)
		plot(mapleRange[c(1,3,47,43),], border=rangeBorder, col=rangeCol, add=T, lwd=rangeLwd)
		plot(grLakes, col="white", add=T)
})
	

