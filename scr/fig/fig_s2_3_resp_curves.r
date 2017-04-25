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


filename = file.path('img', 'figs', 'fig_s2.pdf')
cex.axis = 0.7
cex.xtitle = 1
cex.ytitle = 1
cex.title = 0.85
line.title = 0
xlims = c(-5, 20)
layout.height=c(1,1,.15,1,.15,1,1,1)

dpi = 600
figure.width = 12/2.54
figure.height = 16/2.54
fontsize=10
mar=c(1.5,2,1.5,0.5)
oma = c(2,2,1.5,3)
pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
## png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
## 	file=filename, pointsize=fontsize, res=dpi)
## par(mfrow=c(6, 4), bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
layout(matrix(c(1:8, rep(25,4), 9:12, rep(26,4), 13:24), ncol=4, byrow=T), height=layout.height)
par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)

for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_0_respCurve.rds')))
	if(is.na(info$rc_tymax)) yl = c(0, 1) else yl = c(0, info$rc_tymax)
	plot(rc$temp, rc$col.temp, col=col.color, xlim=xlims, type='l', xlab='', ylab='', ylim=yl, main=plLab, cex.main=cex.title)
## 	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$temp, rev(rc$temp)), c(rc$col.temp.lower, rev(rc$col.temp.upper)), col=col.bgcolor, border=NA)
	lines(rc$temp, rc$ext.temp, col=ext.color) 
	polygon(c(rc$temp, rev(rc$temp)), c(rc$ext.temp.lower, rev(rc$ext.temp.upper)), col=ext.bgcolor, border=NA)
}

mtext("Mean Annual Temperature (°C)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext("Probability of colonization/extinction",col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
## 		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
## 		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
## 		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)

# labels
cex.ti = 1.4
par(xpd=NA)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.84, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.59, "ndc", "user"), "Transitional", 
		srt=270, cex=cex.ti)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.3, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)

dev.off()


filename = file.path('img', 'figs', 'fig_s3.pdf')
## png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
## 	file=filename, pointsize=fontsize, res=dpi)
pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
layout(matrix(c(1:8, rep(25,4), 9:12, rep(26,4), 13:24), ncol=4, byrow=T), height=layout.height)
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_0_respCurve.rds')))
	if(is.na(info$rc_pymax)) yl = c(0, 1) else yl = c(0, info$rc_pymax)
	plot(rc$precip, rc$col.precip, col=col.color, type='l', xlab='', ylab='', ylim=yl, main=plLab, cex.main=cex.title)
## 	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$precip, rev(rc$precip)), c(rc$col.precip.lower, rev(rc$col.precip.upper)), col=col.bgcolor, border=NA)
	lines(rc$precip, rc$ext.precip, col=ext.color) 
	polygon(c(rc$precip, rev(rc$precip)), c(rc$ext.precip.lower, rev(rc$ext.precip.upper)), col=ext.bgcolor, border=NA)
}
mtext("Total Annual Precipitation (mm)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext("Probability of colonization/extinction",col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
## 		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
## 		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
## mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
## 		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)

par(xpd=NA)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.84, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.59, "ndc", "user"), "Transitional", 
		srt=270, cex=cex.ti)
text(grconvertX(0.97, "ndc", "user"), grconvertY(0.3, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)

dev.off()