#!/usr/bin/Rscript
library(coda)
library(sp)
library(raster)
library(rgdal)

## species chosen as the 4 example species
#spNames = c('19481-BET-ALL', '183295-PIC-GLA', '195773-POP-TRE', '19408-QUE-RUB')
spNames = c('18032-ABI-BAL','28731-ACE-SAC')

## spNames = c('18032-ABI-BAL', '28728-ACE-RUB', '28731-ACE-SAC', '19481-BET-ALL', 
## 		'19489-BET-PAP', '19462-FAG-GRA', '32931-FRA-AME', '32929-FRA-PEN', 
## 		'18086-LIR-TUL', '27821-NYS-SYL', '183295-PIC-GLA', '183302-PIC-MAR',  
## 		'183385-PIN-STR', '195773-POP-TRE', '19290-QUE-ALB', 
## 		'19280-QUE-NIG', '19408-QUE-RUB', '19447-QUE-VEL', '19049-ULM-AME')
		

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
	
}



col.colors = colorRampPalette(c('#eff3ff', '#6baed6', '#08519c'), space="rgb")
ext.colors = colorRampPalette(c('#fee5d9', '#fb6a4a', '#a50f15'))
pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
sdm.colors = pres.colors

exp.colors = colorRampPalette(c('#ffffff', '#eff3ff', '#6baed6', '#08519c'))
con.colors = colorRampPalette(c('#ffffff', '#fee5d9', '#fb6a4a', '#a50f15'))
stay.colors = colorRampPalette(c('#ffffff', '#efedf5', '#bcbddc', '#756bb1'))
cat.3colors = c('#1f78b4', '#b2df8a', '#fb9a99')
cat.co.colors = c('#1f78b4', '#fb9a99')
cat.ex.colors = c('#1f78b4', '#b2df8a')
cex.title = 0.5
leg.args = list(cex.axis=0.6)
i = 1
thresh = 0.1

paperwidth = 6.5
dpi = 600
hToWRatio = 1
width = as.integer(dpi*paperwidth)
height = as.integer(width * hToWRatio)
fontsize = 15
png(w=width, h=height, file="img/posterior_maps.png", pointsize=fontsize, res = dpi)
par(mfrow=c(4,3), mar=c(0.5,0,0,0.5), oma=c(0,1,1,0))

for(spName in spNames)
{

	spInfo = spInfoAll[spInfoAll$spName == spName,]
	spLab = bquote(italic(.(spInfo$genus)~.(spInfo$species)))
## 	load(file.path('species', spName, 'res', paste(spName, 'mapRasters.rdata', sep='_')))


	plot(grid.pres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	mtext(spLab, side=2, cex=cex.title)
	if(i == 1) 	mtext("Probability of presence (STM)", cex=cex.title)
	plot(grid.sdm, col=sdm.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("Probability of presence (SDM)", cex=cex.title)
	
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
		
	plot(grid.cat, col=cat.colors, xaxt='n', yaxt='n', legend=FALSE)
	plotbg()
	if(i == 1) 	mtext("Range disequilibrium", cex=cex.title)

	i = i+1
}

dev.off()