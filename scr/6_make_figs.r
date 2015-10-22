#!/usr/bin/env Rscript

## figures produced:
##    img/resp_curves.png (all species)
##    img/posterior_summarys.png (all sp)

library(coda)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

source('scr/stm_functions.r')

# colors
col.color = '#377eb8'
col.bgcolor = '#377eb855'
ext.color = '#e41a1c'
ext.bgcolor = '#e41a1c55'


dpi = 600
figure.width = 6.5
figure.height = 8
filename = file.path('img', 'resp_curves.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(5, 4), bty='n', mar=c(2,2,0,0), mgp=c(1,0.25,0), oma=c(2,2,0,0), tcl=-0.2, cex.axis=0.5)
cex.xtitles = 0.5
cex.ytitle = 0.75
cex.title = 0.5
line.title = -1
rnum = 1
cnum = 1
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_respCurve.rds')))
	plot(rc$temp, rc$col.temp, col=col.color, type='l', xlab='', ylab='', ylim=c(0,info$rc_tymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$temp, rev(rc$temp)), c(rc$col.temp.lower, rev(rc$col.temp.upper)), col=col.bgcolor, border=NA)
	lines(rc$temp, rc$ext.temp, col=ext.color) 
	polygon(c(rc$temp, rev(rc$temp)), c(rc$ext.temp.lower, rev(rc$ext.temp.upper)), col=ext.bgcolor, border=NA)
	if(rnum == 5)
		mtext("Mean Annual Temperature (°)", side = 1, cex=cex.xtitles, line=2)
	plot(rc$precip, rc$col.precip, col=col.color, type='l', xlab='', ylab='', ylim=c(0, info$rc_pymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$precip, rev(rc$precip)), c(rc$col.precip.lower, rev(rc$col.precip.upper)), col=col.bgcolor, border=NA)
	lines(rc$precip, rc$ext.precip, col=ext.color) 
	polygon(c(rc$precip, rev(rc$precip)), c(rc$ext.precip.lower, rev(rc$ext.precip.upper)), col=ext.bgcolor, border=NA)
	if(rnum == 5)
		mtext("Mean Annual Precipitation (mm)", side = 1, cex=cex.xtitles, line=2)
	if(cnum == 2)
	{
		cnum = 1
		rnum = rnum+1
	} else
		cnum = cnum + 1
}
mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
dev.off()





dpi = 600
figure.width = 6.5
figure.height = 15
filename = file.path('img', 'posterior_summaries.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(10, 4), mar=c(1,0,1,0), oma=c(2,1,1,0), tcl=-0.2, cex.axis=0.5)

stmCols = colorRampPalette(c("#ffffff", rev(c("#66c2a4", "#2ca25f", '#006d2c'))), 
		interpolate='spline', space="rgb")
sdmCols = colorRampPalette(c("#ffffff", rev(c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", '#253494'))), 
	interpolate='spline', space="rgb")
stmVarCols = colorRampPalette(c("#ffffff", "#fee5d9", '#fcae91', '#fb6a4a', '#de2d26', 
		'#a50f15'), interpolate='spline', space="rgb")

rdeCols = list(
	cur = 2,
	get_next = function()
	{
		rdeCols$cur <<- rdeCols$cur + 1
		if(rdeCols$cur > length(rdeCols)) rdeCols$cur <<- 3
		rdeCols[[rdeCols$cur]]
	},
	c('#1f78b4', '#b2df8a', '#fb9a99'),
	c('#1f78b4','#a6cee3','#b2df8a'),
	c('#386cb0','#beaed4','#FFE87C')
)

for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]

	spGrid = readRDS(file.path('res','maps',paste0(spName,'_maps.rds')))
	spGrid$stm[spGrid$stm == 0] = NA
	spGrid$sdm[spGrid$sdm == 0] = NA
	plot_sdm(spGrid$stm, spGrid[,1:2], stmCols(100))
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	mtext(plLab, side=2, cex=0.6)
	plot_sdm(spGrid$stm.var, spGrid[,1:2], stmVarCols(100))
	plot_sdm(spGrid$sdm, spGrid[,1:2], sdmCols(100))

	plot_sdm(spGrid$rde, spGrid[,1:2], rdeCols$get_next())
}
dev.off()

