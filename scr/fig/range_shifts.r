library(Hmisc)
library(raster)
library(rgdal)

setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spNames = c('18032-ABI-BAL', '28728-ACE-RUB', '28731-ACE-SAC', '19481-BET-ALL', 
		'19489-BET-PAP', '19462-FAG-GRA', '32931-FRA-AME', '32929-FRA-PEN', 
		'18086-LIR-TUL', '27821-NYS-SYL', '183295-PIC-GLA', '183302-PIC-MAR',  
		'183385-PIN-STR', '195773-POP-TRE', '19290-QUE-ALB', 
		'19280-QUE-NIG', '19408-QUE-RUB', '19447-QUE-VEL', '19049-ULM-AME')

load("dat/map_projections.rdata")

qvals = c(0.25, 0.75)
ns = length(spNames)
stm.mean = sdm.mean = matrix(NA, nrow=ns, ncol=2)
sdm.quant = stm.quant = matrix(NA, ncol=4, nrow=ns)
i = 1
for(spName in spNames)
{
	load(file.path('species', spName, 'res', paste(spName, 'mapRasters.rdata', sep="_")))
	mean.lat.sdm = wtd.mean(coordinates(grid.sdm)[,2], getValues(grid.sdm))
	mean.lon.sdm = wtd.mean(coordinates(grid.sdm)[,1], getValues(grid.sdm))
	mean.lat.stm = wtd.mean(coordinates(grid.pres)[,2], getValues(grid.pres))
	mean.lon.stm = wtd.mean(coordinates(grid.pres)[,1], getValues(grid.pres))
	
	stm.mean[i,] = apply(coordinates(grid.pres), 2, wtd.mean, weights=getValues(grid.pres))
	sdm.mean[i,] = apply(coordinates(grid.sdm), 2, wtd.mean, weights=getValues(grid.sdm))
	stm.quant[i, ] = as.vector(apply(coordinates(grid.pres), 2, wtd.quantile, weights=getValues(grid.pres), probs=qvals))
	sdm.quant[i, ] = as.vector(apply(coordinates(grid.sdm), 2, wtd.quantile, weights=getValues(grid.sdm), probs=qvals))
	
	i = 1 + i
}


xlim=range(c(stm.quant[,1:2], sdm.quant[,1:2])) 
ylim=range(c(stm.quant[,3:4], sdm.quant[,3:4])) 

plot(0,0,type='n',  xlim=xlim, ylim=ylim, xlab="longitude", ylab="latitude")

ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)
lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
# grab specific lakes
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
grLakes = spTransform(grLakes, stmMapProjection)

plotbg = function(txt="", rangeMap=T)
{
	plot(ocean, col="white", add=T)
	plot(grLakes, col="white", add=T)
	mtext(txt, adj=0, cex = 1)
	
	# plot grid lines and labels
	llgridlines(ocean, easts=seq(-90, -40, 10), norths=seq(20,50,10), ndiscr=100, 
			side="ES", offset=1e6, col='#44444444')
	par(xpd=TRUE)
#	westText = c("30°N", "40°N", "50°N")
	eastText = c("20°N", "30°N", "40°N")
#	westCoords = list(rep(-3e5,3), c(0.75e6, 1.9e6, 3e6))
	eastCoords = list(rep(3.5e6,3), c(0.55e6, 1.7e6, 2.9e6))
	southText = c("90°W", "80°W", "70°W")
	southCoords = list(c(0.65e6, 1.7e6, 2.8e6), rep(-0.7e6,3))
	cex.ll = 0.8
	text(eastCoords[[1]], eastCoords[[2]], eastText, cex=cex.ll)
#	text(westCoords[[1]], westCoords[[2]], westText, cex=cex.ll)
	text(southCoords[[1]], southCoords[[2]], southText, cex=cex.ll)
	par(xpd=FALSE)
}

plotbg()

points(stm.mean[,1], stm.mean[,2], pch=16, col='blue',)
points(sdm.mean[,1], sdm.mean[,2], pch=16, col='red')	
## segments(stm.quant[,1], stm.mean[,2], stm.quant[,2],stm.mean[,2], col='blue')
## segments(sdm.quant[,1], sdm.mean[,2], sdm.quant[,2],sdm.mean[,2], col='red')
## segments(stm.mean[,1], stm.quant[,3], stm.mean[,1], stm.quant[,4],, col='blue')
## segments(sdm.mean[,1], sdm.quant[,3], sdm.mean[,1], sdm.quant[,4],, col='red')
arrows(sdm.mean[,1], sdm.mean[,2], stm.mean[,1], stm.mean[,2], length=0.1)
