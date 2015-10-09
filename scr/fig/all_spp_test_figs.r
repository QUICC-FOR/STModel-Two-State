#!/usr/bin/Rscript
library(coda)
library(sp)
library(raster)
library(rgdal)

# get list of all species by reading the directory list
spList = list.files('species')
spInfoAll = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE, colClasses='character')
varNames = readRDS("dat/climVariableNames.rds")

# get base data
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
}


# color ramps



# png settings


spPlot = function(sp)
{
	spInfo = spInfoAll[spInfoAll$spName == sp,]
	spLab = bquote(italic(.(spInfo$genus)~.(spInfo$species)))
	load(file.path('species', sp, 'res', paste(sp, 'mapRasters.rdata', sep='_')))
	e1 = spInfo$env1
	e2 = spInfo$env2
	ylim1 = as.numeric(spInfo$rylim1)
	ylim2 = as.numeric(spInfo$rylim2)
	e1.lab = varNames[which(varNames[,1] == e1),2]
	e2.lab = varNames[which(varNames[,1] == e2),2]

	resp = readRDS(file.path('species', spName, 'res', paste(spName, 'responseCurves.rds', sep='_')))

	# plot colors and settings
	pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
	sdm.colors = pres.colors
	leg.args = list(cex.axis=0.9)
	cex.title = 0.75
	cat.3colors = c('#1f78b4', '#b2df8a', '#fb9a99')
	cat.co.colors = c('#1f78b4', '#fb9a99')
	cat.ex.colors = c('#1f78b4', '#b2df8a')
	thresh = 0.1

	figwidth=18
	figaspect=4.5
	figheight=figwidth/figaspect
	
	dpi = 600
	width = as.integer(dpi*figwidth)
	height = as.integer(dpi*figheight)
	fontsize = 15
	fname = file.path('species', sp, 'img', 'posterior_figs.png')
	png(w=width, h=height, file=fname, pointsize=fontsize, res = dpi)
	
	# make plots
	par(mfrow=c(1,5), mar=c(0.5,0,0,3.5), oma=c(0,1,1,0), bty='o')

	plot(grid.pres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	mtext(spLab, side=2, cex=cex.title)
	mtext("Probability of presence (STM)", cex=cex.title)

	plot(grid.sdm, col=sdm.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	mtext("Probability of presence (SDM)", cex=cex.title)
	
	
	contract = which(getValues(grid.pres) < thresh & getValues(grid.sdm) > thresh)
	expand = which(getValues(grid.sdm) < thresh & getValues(grid.pres) > thresh)
	absent = which(getValues(grid.pres) < thresh & getValues(grid.sdm) < thresh)
	stay = which(!(1:length(grid.pres) %in% c(contract, expand, absent)))
	
	grid.expand = grid.pres - grid.sdm
	grid.expand[c(contract, stay, absent)] = NA
	grid.contract =  grid.sdm - grid.pres
	grid.contract[c(expand, stay, absent)] = NA
	grid.stay = grid.pres - grid.sdm
	grid.stay[c(expand, contract, absent)] = NA
	
	grid.cat = 0 * grid.pres
	grid.cat[expand] = 1
	grid.cat[contract] = 2
	grid.cat[absent] = NA
	
	if(length(contract) == 0) {
		cat.colors = cat.ex.colors
	} else if(length(expand) == 0) {
		cat.colors = cat.co.colors
	} else {
		cat.colors = cat.3colors
	}
		
	par(mar=c(0.5,0,0,0.5))
	plot(grid.cat, col=cat.colors, xaxt='n', yaxt='n', legend=FALSE)
	plotbg()
	mtext("Range disequilibrium", cex=cex.title)

	par(mar=c(3.5,3.5,0,0.9), bty='n', mgp=c(2, 0.5, 0))
	with(resp,
	{
		plot(env1, col1.mean, xlab=e1.lab, ylab="Probability", ylim=c(0,ylim1), col='blue', type='l')
		polygon(c(env1, rev(env1)), c(col1.lo, rev(col1.up)), col="#0000FF22", border=NA)
		lines(env1, ext1.mean, col='red')
		polygon(c(env1, rev(env1)), c(ext1.lo, rev(ext1.up)), col="#FF000022", border=NA)

		plot(env2, col2.mean, xlab=e2.lab, ylab="", ylim=c(0,ylim2), col='blue', type='l')
		polygon(c(env2, rev(env2)), c(col2.lo, rev(col2.up)), col="#0000FF22", border=NA)
		lines(env2, ext2.mean, col='red')
		polygon(c(env2, rev(env2)), c(ext2.lo, rev(ext2.up)), col="#FF000022", border=NA)
		xpd = par()$xpd
		par(xpd=NA)
		par(xpd=xpd)
	})
	dev.off()
}

for(spName in spList)
{
	# try to read map data & make plots
	try(spPlot(spName))


}


