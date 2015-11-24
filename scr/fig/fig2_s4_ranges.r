#!/usr/bin/env Rscript

## figures produced:
##    img/figs/fig2.png

source('scr/stm_functions.r')
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')
suppressWarnings(dir.create(file.path('img', 'figs'), recursive=TRUE))

arg = commandArgs(TRUE)
mod = if(length(arg) > 0) arg[1] else '0'

dpi = 300
figure.width = 12
figure.height = 9.6
filename = file.path('img', 'figs', 'fig2')
fontsize=7

mainTxtSps = c("28731-ACE-SAC", "19462-FAG-GRA", "32931-FRA-AME", "19408-QUE-RUB",
		"19481-BET-ALL", "183375-PIN-RES", "183385-PIN-STR", "22463-POP-GRA",
		"18032-ABI-BAL",  "19489-BET-PAP", "183295-PIC-GLA", "195773-POP-TRE")

## png(width=figure.width, height=figure.height, units='cm',
## 	file=paste0(filename, '.png'), pointsize=fontsize, res=dpi)
pdf(width=figure.width/2.54, height=figure.height/2.54, file=paste0(filename,'.pdf'), pointsize=fontsize)
par(mfrow=c(3, 4), mar=c(0,0,2,0), oma=c(1,1,0,2), tcl=-0.2, cex.axis=0.5)
## 
rdeCols = c('#1f78b4', '#b2df8a', '#fb9a99')
## 
for(spName in mainTxtSps)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))

	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	plot.sdm(spGrid[,c('x', 'y', 'rde')], sdm.col=rdeCols, box=TRUE, axes=FALSE, main=plLab)
}

# legend
labs = c("Continued presence", "Colonization credit", "Extinction debt")
## plot(0,0,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), type='n')

legend(grconvertX(0.05, "npc", "user"), grconvertY(0.35, "npc", "user"), fill=rdeCols, 
		legend=labs, cex=1, bg='white')

# labels
cex.ti = 1.5
par(xpd=NA)
text(grconvertX(0.975, "ndc", "user"), grconvertY(0.8, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(0.975, "ndc", "user"), grconvertY(0.5, "ndc", "user"), "Transitional", 
		srt=270, cex=cex.ti)
text(grconvertX(0.975, "ndc", "user"), grconvertY(0.2, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)

dev.off()


siSps = mainTxtSps = c("28728-ACE-RUB", "32945-FRA-NIG", "19287-QUE-MAC", "183397-TSU-CAN",
		"183412-LAR-LAR", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN", "505490-THU-OCC")

figure.width = 12
figure.height = 6
pdf(width=figure.width/2.54, height=figure.height/2.54, file='img/figs/fig_s4.pdf', pointsize=fontsize)
par(mfrow=c(2, 5), mar=c(0,0,2,0), oma=c(1,1,0,2), tcl=-0.2, cex.axis=0.5)
for(spName in siSps)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))

	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	plot.sdm(spGrid[,c('x', 'y', 'rde')], sdm.col=rdeCols, box=FALSE, axes=TRUE, main=plLab)
	if(spName == "183397-TSU-CAN")
	{
		plot(0,0, xaxt='n', yaxt='n', bty='n', xlab='', ylab='', type='n')
		labs = c("Continued presence", "Colonization credit", "Extinction debt")
		legend(grconvertX(0.05, "npc", "user"), grconvertY(0.5, "npc", "user"), fill=rdeCols, 
				legend=labs, cex=0.8, bg='white')
	}
}

text(grconvertX(0.975, "ndc", "user"), grconvertY(0.8, "ndc", "user"), "Temperate", 
		srt=270, cex=cex.ti)
text(grconvertX(0.975, "ndc", "user"), grconvertY(0.2, "ndc", "user"), "Boreal", 
		srt=270, cex=cex.ti)
dev.off()
