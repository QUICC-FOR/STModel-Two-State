#!/usr/bin/env Rscript
# make a figure comparing posterior parameters from INLA and mcmc

library("INLA")
source("scr/inla/inla_functions.r")
speciesInfo = read.csv('dat/speciesInfo.csv')

spAll <- readRDS("dat/speciesList.rds")
spList <- commandArgs()
spList <- spList[spList %in% spAll]
if(length(spList) == 0) spList <- spAll

#####
spName <- spList[1]

pSummary <- function(pars) cbind(mean=colMeans(pars), lower=apply(pars, 2, quantile, 0.05), 
																 upper=apply(pars, 2, quantile, 0.95))
sep <- 0.13
xPos <- c((1:5) - sep, 1:5, (1:5)+sep)
pch <- 20
catCols <- c('#0272fd', '#d78c39', '#aa004b')
cols <- rep(catCols, each=5)
legendLabs <- c("Full Model", "5-year nonspatial", "5-year spatial")

figure.width = 12/2.54
figure.height = 16/2.54
filename = file.path('img', 'inla', 'inla_params.pdf')
fontsize=10
layout.height=c(1,1,.15,1,.15,1,1,1)
mar=c(1, 0, 2, 2)
oma = c(2.5, 3, 0, 1)
cex.axis = 0.7
cex.xtitle = 1
cex.ytitle = 1
cex.title = 0.8


pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
layout(matrix(c(1:8, rep(25,4), 9:12, rep(26,4), 13:24), ncol=4, byrow=T), height=layout.height)
par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
for(spName in spList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))
	resDir <- file.path('res', 'inla', spName)
	inla_samps <- readRDS(file.path(resDir, "inla_samples.rds"))
	mcmc <- read.csv(file.path('res', '/mcmc/', spName, '0', '/posterior.csv'))
	mcmc <- mcmc[seq(1, nrow(mcmc), length.out=1000),8:12]  ## extinction parameters only
	mcPars <- pSummary(mcmc)
	inPars <- lapply(inla_samps$samples, function(x) pSummary(t(x[1:5,])))
	params <-do.call(rbind, c(list(mcPars), inPars))
	plot(xPos, params[,1], pch=pch, xlim=c(0.5,5.5), ylim=range(params),
			 col=cols, xaxt='n', xlab="", ylab="")
	title(plLab, cex.main=cex.title, line=0.5)
	segments(xPos, params[,2], xPos, params[,3], col=cols)
	axis(side=1, at=1:5, labels=c('i', 'T', 'P', expression(T^2), expression(P^2)))
}
plot.new()
legend("topleft", col=catCols, lwd=1.5, legend=legendLabs, bty='n', cex=0.8)
mtext("Parameter", side=1, outer=TRUE, cex=cex.xtitle, line=1)
mtext("Value", side=2, outer=TRUE, cex=cex.ytitle, line=2)
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
