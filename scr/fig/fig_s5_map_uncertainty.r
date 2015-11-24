#!/usr/bin/env Rscript

## figures produced:
##    img/figs/fig_s4.png


library(RColorBrewer)
varCol = colorRampPalette(c('#ffffff', brewer.pal(7, "YlOrRd")), interpolate='spline', bias=1)(200)

source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

stmColors = colorRampPalette(brewer.pal(5, "YlGnBu"),
		interpolate='spline', bias=1, space="rgb")(200)

leg.title.cex = 1
leg.axis.cex = 1
cex.ti = 1.5
temperate.ypos = 0.82
transitional.ypos = 0.58
boreal.ypos = 0.35
title.xpos = 0.96
	
figure.width = 16
figure.height = 24
fname = file.path('img', 'figs', 'fig_s5.png')


fname = file.path('img', 'figs', paste0('fig_s5.png'))
## pdf(file=paste0(fname, '.pdf'), width=figure.width/2.54, height=figure.height/2.54, pointsize=9)
png(file=fname, width=figure.width, height=figure.height, pointsize=9, res=300, units='cm')
## pdf(file=fname, width=figure.width/2.54, height=figure.height/2.54, pointsize=9)
par(mfrow=c(6, 4), mar=c(0,0,2,0), oma=c(1,1,0,2), tcl=-0.2)
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))

	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	plot.sdm(with(spGrid, data.frame(x, y, sqrt(stm*(1-stm)))), varCol, legend=FALSE, axes=FALSE, box=TRUE, main=plLab)
## 	mtext(plLab, side=3, cex=0.6)
}
plot.new()
plot.sdm(with(spGrid, data.frame(x, y, sqrt(stm*(1-stm)))), varCol, legend=TRUE, legend.only=TRUE, 
		plot.ocean=FALSE, legend.width=6, legend.shrink=1, horizontal=TRUE, 
		legend.args=list(text='Posterior standard deviation', side=3, line=0.5, 
		cex=leg.title.cex, cex.axis = leg.axis.cex), smallplot=c(0.1, 0.9, 0.2, 0.35), legend.cex=1)


plot.new()
par(xpd=NA)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(temperate.ypos, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(transitional.ypos, "ndc", "user"), "Transitional", 
		srt=270, cex=cex.ti)
text(grconvertX(title.xpos, "ndc", "user"), grconvertY(boreal.ypos, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)
dev.off()
