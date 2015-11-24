#!/usr/bin/env Rscript

## figures produced:
##    img/figs/fig2.png

library(coda)

source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
suppressWarnings(dir.create(file.path('img', 'figs'), recursive=TRUE))
# speciesInfo = read.csv('dat/speciesInfo.csv')

arg = commandArgs(TRUE)
mod = if(length(arg) > 0) arg[1] else '0'

dpi = 600
figure.width = 6.5
figure.height = 4.5
filename = file.path('img', 'figs', paste0('fig2_', mod, '.png'))
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(3, 4), mar=c(0,0,1,0), oma=c(1,1,0,0), tcl=-0.2, cex.axis=0.5)

rdeCols = c('#1f78b4', '#b2df8a', '#fb9a99')

for(spName in speciesList)
{

	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))

	spGrid = readRDS(file.path('res','maps',paste0(spName,'_', mod, '_maps.rds')))
	
	# get the calibration range
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
	latLim = range(calibDat$lat)
	lonLim = range(calibDat$lon)
	
## 	# rde is fucked up for some reason; fix it here temporarily
## 	spGrid = within(spGrid,
## 	{
## 		rde[rde.present >= rde.contract & rde.present >= rde.expand ] = 0
## 		rde[rde.expand >= rde.contract & rde.expand > rde.present ] = 1
## 		rde[rde.contract > rde.present & rde.contract > rde.expand ] = 2
## 		rde[sdm < 0.1 & stm < 0.1 ] = NA
## 		
## 		# restrict to calibration range
## 		rde[lon < lonLim[1] | lon > lonLim[2] | lat < latLim[1] | lat > latLim[2]] = NA
## 		if(spName == "183302-PIC-MAR")
## 			rde[lat < (latLim[1]+1.9)] = NA
## 	})

	plot_sdm(spGrid$rde, spGrid[,1:2], sdm.col=rdeCols)
	mtext(plLab, side=3, cex=0.6)
}


# legend
par(xpd=NA)

labs = c("Continued presence", "Colonization lag", "Extinction debt")
plot(0,0,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), type='n')
legend(0.1, 0.9, fill=rdeCols, legend=labs)

dev.off()

