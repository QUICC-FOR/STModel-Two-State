#!/usr/bin/env Rscript
# compute the posterior distribution of the response curves for inla models
library("INLA")
source("scr/inla/inla_functions.r")
speciesInfo = read.csv('dat/speciesInfo.csv')

spAll <- readRDS("dat/speciesList.rds")
spList <- commandArgs()
spList <- spList[spList %in% spAll]
if(length(spList) == 0) spList <- spAll

climScale = readRDS('dat/clim/climate_scaling.rds')
climDat = readRDS('dat/clim/plotClimate_scaled.rds')

rcRes = 500

# set ranges and non-constant x-values for the response curves
trange <- c(-5, 25)
prange <- c(500, 2000)
tempUnscaled <- seq(trange[1], trange[2], length.out=rcRes)
precipUnscaled <- seq(prange[1], prange[2], length.out=rcRes)
rc.env1 <- scale(tempUnscaled, center = climScale$center['annual_mean_temp'], 
								 scale = climScale$scale['annual_mean_temp'])
rc.env2 <- scale(precipUnscaled, center = climScale$center['tot_annual_pp'], 
								 scale = climScale$scale['tot_annual_pp'])

cols <- list(mcmc='#0272fd', nonspatial='#d78c39', spatial='#aa004b')
legendLabs <- c("Full Model", "5-year nonspatial", "5-year spatial")
filename = file.path('img', 'inla', 'inla_rc_temp.pdf')
cex.axis = 0.7
cex.xtitle = 1
cex.ytitle = 1
cex.title = 0.85
line.title = 0
xlims = c(-5, 20)
layout.height=c(1,1,.15,1,.15,1,1,1)
figure.width = 12/2.54
figure.height = 16/2.54
fontsize=10
lwd=1
mar=c(1.5,2,1.5,0.5)
oma = c(2,2,1.5,3)
pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
layout(matrix(c(1:8, rep(25,4), 9:12, rep(26,4), 13:24), ncol=4, byrow=T), height=layout.height)
par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
for(spName in spList)
{
	resDir <- file.path('res', 'inla', spName)
	dir.create(resDir, showWarnings = FALSE)
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))

	# get rc constants where the species is present
	presDat = readRDS(file.path('dat', 'presence', paste0(spName,'_presence.rds')))
	presDat = merge(presDat, climDat)
  rc.env2.c = info$rc_pval
  if(is.na(rc.env2.c))
	  rc.env2.c = mean(presDat[presDat[,spName] == 1,'tot_annual_pp'])

	inla_samps <- readRDS(file.path(resDir, "inla_samples.rds"))
	mcmc <- read.csv(file.path('res', '/mcmc/', spName, '0', '/posterior.csv'))
	mcmc.col <- mcmc[seq(1, nrow(mcmc), length.out=1000),1:5]  ## colonization parameters only
	mcmc <- mcmc[seq(1, nrow(mcmc), length.out=1000),8:12]  ## extinction parameters only
	lp.temp <- inla_lp(inla_samps$samples, mcmc, NULL, rc.env1, rc.env2.c)
	lp.col.temp <- cbind(m=1, t=rc.env1, p=rc.env2.c, t2=rc.env1^2, p2=rc.env2.c^2) %*% t(as.matrix(mcmc.col))
	rate.temp <- lapply(lp.temp, plogis)
	rate.col.temp <- plogis(lp.col.temp)
	rc.temp <- lapply(rate.temp, function(x) {
		data.frame(x = tempUnscaled,
							 y = rowMeans(x),
							 lower = apply(x, 1, quantile, 0.05),
							 upper = apply(x, 1, quantile, 0.95))})
	lims <- inla_rc_lims(spName)
	plot(tempUnscaled, rc.temp[[1]]$y, type='n', xlim=lims$x, ylim=lims$y, main=plLab, cex.main=cex.title,
			 xlab="", ylab="")
	mapply(function(rc, col) {
		lines(rc$x, rc$y, col=col, lwd=lwd)
		polygon(c(rc$x, rev(rc$x)), c(rc$lower, rev(rc$upper)), border=NA, col=paste0(col, '33'))
	}, rc.temp, cols)
	lines(rc.temp[[1]]$x, rowMeans(rate.col.temp), col='#666666', lwd=lwd)
}
plot.new()
legend("topleft", col=unlist(cols), lwd=1.5, legend=legendLabs, bty='n', cex=0.8)
mtext("Mean Annual Temperature (°C)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext("Probability of colonization/extinction",col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)

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







filename = file.path('img', 'inla', 'inla_rc_precip.pdf')
pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
layout(matrix(c(1:8, rep(25,4), 9:12, rep(26,4), 13:24), ncol=4, byrow=T), height=layout.height)
par(bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
for(spName in spList)
{
	resDir <- file.path('res', 'inla', spName)
	dir.create(resDir, showWarnings = FALSE)
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1), '.'))~.(as.character(info$species))))
	
	# get rc constants where the species is present
	presDat = readRDS(file.path('dat', 'presence', paste0(spName,'_presence.rds')))
	presDat = merge(presDat, climDat)
	rc.env1.c = info$rc_tval
	if(is.na(rc.env1.c))
		rc.env1.c = mean(presDat[presDat[,spName] == 1,'annual_mean_temp'])
	
	inla_samps <- readRDS(file.path(resDir, "inla_samples.rds"))
	mcmc <- read.csv(file.path('res', '/mcmc/', spName, '0', '/posterior.csv'))
	mcmc.col <- mcmc[seq(1, nrow(mcmc), length.out=1000),1:5]  ## colonization parameters only
	mcmc <- mcmc[seq(1, nrow(mcmc), length.out=1000),8:12]  ## extinction parameters only
	lp.precip <- inla_lp(inla_samps$samples, mcmc, NULL, rc.env1.c, rc.env2)
	lp.col.precip <- cbind(m=1, t=rc.env1.c, p=rc.env2, t2=rc.env1.c^2, p2=rc.env2^2) %*% t(as.matrix(mcmc.col))
	rate.precip <- lapply(lp.precip, plogis)
	rate.col.precip <- plogis(lp.col.precip)
	rc.precip <- lapply(rate.precip, function(x) {
		data.frame(x = precipUnscaled,
							 y = rowMeans(x),
							 lower = apply(x, 1, quantile, 0.05),
							 upper = apply(x, 1, quantile, 0.95))})
	lims <- inla_rc_lims(spName)
	plot(precipUnscaled, rc.precip[[1]]$y, type='n', xlim=lims$xp, ylim=lims$yp, main=plLab, cex.main=cex.title,
			 xlab="", ylab="")
	mapply(function(rc, col) {
		lines(rc$x, rc$y, col=col, lwd=lwd)
		polygon(c(rc$x, rev(rc$x)), c(rc$lower, rev(rc$upper)), border=NA, col=paste0(col, '33'))
	}, rc.precip, cols)
	lines(rc.precip[[1]]$x, rowMeans(rate.col.precip), col='#666666', lwd=lwd)
}
plot.new()
legend("topleft", col=unlist(cols), lwd=1.5, legend=legendLabs, bty='n', cex=0.8)
mtext("Total Annual Precipitation (mm)", side=1, outer=TRUE, cex=cex.xtitle, line=0)
mtext("Probability of colonization/extinction",col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)

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

