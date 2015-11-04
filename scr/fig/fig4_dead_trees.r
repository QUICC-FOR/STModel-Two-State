#!/usr/bin/env Rscript
library(rstan)
rdeCols = c('#1f78b4', '#fb9a99')
speciesInfo = read.csv('dat/speciesInfo.csv')
stnResultPath = file.path('res', 'deadStanMod.rds')

stnResults = readRDS(stnResultPath)

samples = extract(stnResults$results)
predictions = with(samples, {
	list(inrange=plogis(alpha), outrange=plogis(alpha+beta))
})
	
plotPoints = do.call(rbind, lapply(names(predictions), function(j)
{
	do.call(rbind, lapply(1:ncol(predictions[[j]]), function(i) {
		data.frame(
			species=stnResults$species$name[stnResults$species$value == i],
			mean = mean(predictions[[j]][,i]),
			sd = sd(predictions[[j]][,i]),
			upper = quantile(predictions[[j]][,i], 0.975),
			lower = quantile(predictions[[j]][,i], 0.025),
			type = j	
	)}))
}))

cex.axis=0.75
gap=0.3
boxwidth = 0.4
between=5
cex.lab=0.4   # species labels
cex.label = 0.7 # axis labels
srt=0
legend.cex=0.5
nspecies = ncol(predictions[[1]])
ylabpos = rep(0, nspecies)
cols=paste0(rep(rdeCols, each=nspecies), 'aa')
xlocs = rep(rep(seq(0, by=between, length.out=nspecies/2), 2) + rep(c(0,gap), each=nspecies/2),2)
xlim = range(xlocs) + 0.5*c(-between, between)
spInd = list(borealSp = c(1,2,6,7,8), tempSp = c(3,4,5,9,10))
spInd = lapply(spInd, function(x) c(x, x+10))
plotPoints$xloc = NA
plotPoints$xloc[c(spInd$borealSp, spInd$tempSp)] = xlocs
## xlocs = rep(seq(0, by=between, length.out=nspecies), 2) + rep(c(0,gap), each=nspecies)
pdf(width=3.5, height=6, file=file.path('img', 'figs', 'fig4.pdf'))
par(mfrow=c(2,1), mar=c(0.5,2.6,0.5,0.5), cex.axis=cex.axis, mgp=c(1.3,0.4,0), cex.lab=cex.label, tcl=-0.2)
lapply(spInd, function(inds)
{
	plot(0,0, type='n', xlim=xlim, ylim=c(0,0.25), xaxt='n', ylab="Proportion dead", xlab="")

	for(i in inds)
	{
		info = speciesInfo[speciesInfo$spName == as.character(plotPoints$species[i]),]
		plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1),'.'))~.(as.character(info$species))))

		xl = c(plotPoints$xloc[i] - boxwidth, plotPoints$xloc[i]+boxwidth)
		yl = c(plotPoints$lower[i], plotPoints$upper[i])
		polygon( c(xl, rev(xl)), rep(yl, each=2), col=cols[i], border=NA)
		if(i <= nspecies)
			text(plotPoints$xloc[i] + gap/2, ylabpos[i], plLab, cex=cex.lab, srt=srt)
	}
	points(plotPoints$xloc[inds], plotPoints$mean[inds], col='black', pch=16, cex=0.5)
	if(1 %in% inds)
	{
		par(font=2)
		legend(xlocs[length(xlocs)-1], 0.25, legend=c("In Range", 
			"Extinction Debt"), fill=rdeCols, cex=legend.cex, text.col=rdeCols, bty='n')
	}
})
dev.off()