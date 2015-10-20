library(coda)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

source('scr/stm_functions.r')


## response curves can go on all one fig
## do all species on one plot, draw a box around each pair of plots
## use this to accomplish the boxing:

## layout(matrix(c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16), 4, 4, byrow=TRUE))
## replicate(16, hist(rnorm(100)))
## par(xpd=NA)
## rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
##      grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
## rect( grconvertX(0.005, from='ndc'), grconvertY(0.005, from='ndc'),
##      grconvertX(0.495, from='ndc'), grconvertY(0.495, from='ndc'))
## rect( grconvertX(0.505, from='ndc'), grconvertY(0.505, from='ndc'),
##      grconvertX(0.995, from='ndc'), grconvertY(0.995, from='ndc'))
## rect( grconvertX(0.505, from='ndc'), grconvertY(0.005, from='ndc'),
##      grconvertX(0.995, from='ndc'), grconvertY(0.495, from='ndc'))







for(spName in speciesList)
{

	## plots here; dragons!


	paperwidth = 10
	dpi = 600
	hToWRatio = 0.4
	width = as.integer(dpi*paperwidth)
	height = as.integer(width * hToWRatio)
	fontsize = 15
	stPres = make_raster(spGrid$stm, spGrid[,1:2], P4S.latlon, stmMapProjection)
	png(w=width, h=height, file=paste0("img/", spName, "_posterior_maps.png"), pointsize=fontsize, res = dpi)
	par(mfrow=c(1,3))
	pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
	plot(stPres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1.00000001))
	plotbg()
	
	sdPres = make_raster(spGrid$sdm, spGrid[,1:2], P4S.latlon, stmMapProjection)
	plot(sdPres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), legend=FALSE)
	plotbg()
	
	
	cat.colors = c('#1f78b4', '#b2df8a', '#fb9a99')
	diseq = make_raster(spGrid$diseq, spGrid[,1:2], P4S.latlon, stmMapProjection)
	plot(diseq, col=cat.colors, xaxt='n', yaxt='n', legend=FALSE)
	plotbg()

dev.off()
		
		
	
}

