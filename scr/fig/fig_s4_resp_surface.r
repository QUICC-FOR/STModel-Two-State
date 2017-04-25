#!/usr/bin/env Rscript

## figures produced:
##    img/figs/temp.png
##    img/figs/precip.png

## library(coda)
library(fields)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

## source('scr/stm_functions.r')

respSurface <- readRDS("res/resp_surface.rds")
climScale = readRDS('dat/clim/climate_scaling.rds')


# colors

filename = file.path('img', 'figs', 'fig_s4.png')
cex.axis = 0.6
cex.xtitle = 0.85
cex.ytitle = 0.85
cex.title = 0.85
xlims = c(-5, 20)
layout.height=c(1,1,.15,1,.15,1,1,1)

xaxt = c(rep('n', 16), rep('s', 5))
names(xaxt) = speciesList
yaxt = c(rep(c('s', rep('n', 3)), 5), 's')
names(yaxt) = speciesList

dpi = 600
figure.width = 12/2.54
figure.height = 16/2.54
fontsize=10
lwd=0.4
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(6, 4), bty='n', mar=c(1.5,0.5,1,3), oma=c(2,4,1,1.5), mgp=c(1,0.25,0), tcl=-0.2, cex.axis=cex.axis, las=1, cex.main=cex.title)
## par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
ncols = 100
cols <- colorRampPalette(c('#C70608', '#ffffff', '#217FCC'))(ncols)
## cols <- colorRampPalette(c('#ca0020','#f4a582','#ffffff','#92c5de','#0571b0'))(ncols)

invisible(lapply(speciesList, function(spName) {
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))
	if(is.na(info$rc_tymax)) yl = c(0, 1) else yl = c(0, info$rc_tymax)

	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'_stm_calib.rds')))
	temp = respSurface$y
	precip = respSurface$x
	z = respSurface$z[[spName]]
	if(spName == "183302-PIC-MAR") {
		mask <- (z > 0) & matrix(x > 8, nrow=nrow(z), ncol=ncol(z))
		z[mask] <- -0.001
	}
		
	breaks = c(seq(min(z), 0, length.out=ncols/2 + 1), seq(0, max(z), length.out=ncols/2 + 1))
	breaks <- breaks[-(ncols/2 + 1)]
	# x and y are switched in image.plot :-/
	image.plot(x=temp, y=precip, z=z, main=plLab, breaks=breaks, nlevel=ncols, 
			col=cols, xlab='', ylab='', xlim=c(-5, 19), xaxt=xaxt[spName], yaxt=yaxt[spName])
	## 	xx <- (calibDat$tot_annual_pp * climScale$scale['tot_annual_pp']) + climScale$center['tot_annual_pp']
	## 	yy <- (calibDat$annual_mean_temp * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp']
	## points(xx, yy, pch='.', cex=0.7, col='#66666666')
	contour(x=temp, y=precip, z=z, add=TRUE, levels=0, drawlabels=FALSE, lwd=lwd)
}))
mtext("Mean Annual Temperature (°C)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext("Total Annual Precipitation (mm)", side=2, outer=TRUE, cex=cex.ytitle, line=2, las=0)

dev.off()





