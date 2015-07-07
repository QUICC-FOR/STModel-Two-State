#!/usr/bin/Rscript
library(coda)
library(sp)
library(raster)
library(rgdal)

spNames = c('18032-ABI-BAL', '28728-ACE-RUB', '28731-ACE-SAC', '19462-FAG-GRA', 
		'32931-FRA-AME', '183295-PIC-GLA', '183385-PIN-STR', '195773-POP-TRE', 
		'19290-QUE-ALB', '19447-QUE-VEL')
spInfoAll = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE, colClasses='character')

load('dat/map_projections.rdata')

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







col.colors = colorRampPalette(c('#eff3ff', '#6baed6', '#08519c'), space="rgb")
ext.colors = colorRampPalette(c('#fee5d9', '#fb6a4a', '#a50f15'))
pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
sdm.colors = pres.colors
#sdm.colors = colorRampPalette(c("#ffffff" "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
#			interpolate='spline', space="rgb")


pdf(w=6.5, h=8.2, file="img/posterior_maps.pdf")
par(mfrow=c(5,4), mar=c(2,2,2,1), oma=c(0,0,0,2))
cex.title = 0.5
leg.args = list(cex.axis=0.6)
i = 1
for(spName in spNames)
{

	spInfo = spInfoAll[spInfoAll$spName == spName,]
	spLab = bquote(italic(.(spInfo$genus)~.(spInfo$species)))
	load(file.path('species', spName, 'res', paste(spName, 'mapRasters.rdata', sep='_')))


	plot(grid.c, col=col.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) mtext("colonization probability", cex=cex.title)
	mtext(spLab, side=2, cex=cex.title)
	plot(grid.e, col=ext.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("extinction probability", cex=cex.title)
	plot(grid.pres, col=sdm.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("Probability of presence (CE Model)", cex=cex.title)
	plot(grid.sdm, col=sdm.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("Probability of presence (SDM)", cex=cex.title)
	
	i = i+1

}
dev.off()

