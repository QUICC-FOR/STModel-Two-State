#!/usr/bin/env Rscript

## figures produced:
##    img/resp_curves.png (all species)
##    img/posterior_summarys.png (all sp)

library(coda)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
models=c('0', 'g')
source('scr/stm_functions.r')

# colors
col.color = '#377eb8'
col.bgcolor = '#377eb855'
ext.color = '#e41a1c'
ext.bgcolor = '#e41a1c55'


dpi = 600
figure.width = 6.5
figure.height = 8
fontsize=12

for(mod in models)
{
	filename = paste0('img/resp_curves_', mod, '.png')
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
		rc = readRDS(file.path('res','resp_curve', paste0(spName, '_', mod, '_respCurve.rds')))
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
}