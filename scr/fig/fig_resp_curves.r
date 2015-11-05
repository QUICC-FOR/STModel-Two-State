#!/usr/bin/env Rscript

## figures produced:
##    img/figs/temp.png
##    img/figs/precip.png

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
figure.width = 5
figure.height = 7
filename = file.path('img', 'figs', 'temp.png')
fontsize=12
mar=c(2,2,0.5,0.5)
cex.axis = 0.7
cex.xtitle = 0.75
cex.ytitle = 0.75
cex.title = 0.5
line.title = -1
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(5, 2), bty='n', mar=mar, mgp=c(1,0.25,0), oma=c(2,2,0,0), tcl=-0.2, cex.axis=cex.axis)
xlims = list()
xlims[["28728-ACE-RUB"]] = c(-5,20)
xlims[["28731-ACE-SAC"]] = c(0,15)
xlims[["19290-QUE-ALB"]] = c(0,20)
xlims[["19408-QUE-RUB"]] = c(0,20)
xlims[["19049-ULM-AME"]] = c(0,25)
xlims[["18032-ABI-BAL"]] = c(-5,10)
xlims[["19489-BET-PAP"]] = c(-5,10)
xlims[["195773-POP-TRE"]] = c(-5,10)
xlims[["19481-BET-ALL"]] = c(-5,15)
xlims[["183302-PIC-MAR"]] = c(-5,10)

for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_respCurve.rds')))
	plot(rc$temp, rc$col.temp, col=col.color, xlim=xlims[[spName]], type='l', xlab='', ylab='', ylim=c(0,info$rc_tymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$temp, rev(rc$temp)), c(rc$col.temp.lower, rev(rc$col.temp.upper)), col=col.bgcolor, border=NA)
	lines(rc$temp, rc$ext.temp, col=ext.color) 
	polygon(c(rc$temp, rev(rc$temp)), c(rc$ext.temp.lower, rev(rc$ext.temp.upper)), col=ext.bgcolor, border=NA)
}

mtext("Mean Annual Temperature (°C)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
dev.off()


filename = file.path('img', 'figs', 'precip.png')
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(5, 2), bty='n', mar=c(2,2,0,0), mgp=c(1,0.25,0), oma=c(2,2,0,0), tcl=-0.2, cex.axis=cex.axis)
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_respCurve.rds')))
	plot(rc$precip, rc$col.precip, col=col.color, type='l', xlab='', ylab='', ylim=c(0, info$rc_pymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$precip, rev(rc$precip)), c(rc$col.precip.lower, rev(rc$col.precip.upper)), col=col.bgcolor, border=NA)
	lines(rc$precip, rc$ext.precip, col=ext.color) 
	polygon(c(rc$precip, rev(rc$precip)), c(rc$ext.precip.lower, rev(rc$ext.precip.upper)), col=ext.bgcolor, border=NA)
}
mtext("Mean Annual Precipitation (mm)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
dev.off()