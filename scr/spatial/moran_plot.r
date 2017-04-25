#!/usr/bin/env Rscript

## produce a plot of moran's I for all species

# needed files and directories
# res/moran/*

if(interactive()) {
	setwd("/Users/mtalluto/Documents/work/projects_active/stm/STModel-Two-State")
} else {
	arg = commandArgs(TRUE)
	spName = arg[1]
}
speciesInfo = read.csv('dat/speciesInfo.csv')


cex.axis = 0.7
mar=c(2,2,0,0)
oma = c(2,2,0.5,0.5)
figure.width = 16/2.54
figure.height = 20/2.54
mar=c(1.5,2,1.5,0.5)
cex.lab=0.9
fontsize=10

spList = readRDS("dat/speciesList.rds")
pdf(width=figure.width, height=figure.height, file="img/moran.pdf", pointsize=fontsize)
## pdf(file="img/moran.pdf", w=figure.width/2.54, h=figure.width/2.54)
par(mfrow=c(6,4), bty='n', mar=mar, mgp=c(1,0.25,0), oma=oma, tcl=-0.2, cex.axis=cex.axis)
## par(mfrow=c(6,4), mar=c(1,1,1,1), oma=c(1.5,3,0,0), mgp=c(1,0.5,0), tcl=-0.2)
sapply(spList, function(spName) {
	tryCatch({
		info = speciesInfo[speciesInfo$spName == spName,]
		plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
		moranMean <- readRDS(paste0('res/moran/', spName, '_moranMean.rds'))
		moranSamp <- readRDS(paste0('res/moran/', spName, '_moranSamp.rds'))
		plot(moranMean[,4], moranMean[,1], pch=20, xlab="", ylab="", ylim=c(-0.1,0.7), type='b', lty=2)
		abline(h=0, lty=2)
		segments(moranSamp[,1], apply(moranSamp[,-1], 1, quantile, 0.95), moranSamp[,1], 
				 apply(moranSamp[,-1], 1, quantile, 0.05))
		text(max(moranMean[,4]), 0.6, plLab, pos=2, cex=cex.lab)
	}, error=function(e) warning(e))
})
mtext("Lag Distance (km)", side=1, outer=TRUE, line=0)
mtext("Moran's I", side=2, outer=TRUE, line=0.5)
dev.off()









## 
## source("scr/stm_functions.r")
## lags <- 1000 * seq(0, 700, 25)
## numCores = 16
## parSampleSize <- 1000
## datSampleSize <- 20000
## numIters <- 100
## 
## # get data
## datAll <- readRDS(paste0("dat/stm_calib/", spName, "_stm_calib.rds"))
## datAll$interval = datAll$year2 - datAll$year1
## 
## # get par samples and subsample
## samples <- readRDS(paste0("res/posterior/", spName, "_0_samples.rds"))
## samples <- samples[sample(nrow(samples), parSampleSize),]
## 
## # subsample data
## dat <- datAll[sample(nrow(datAll), datSampleSize),]
## 
## resids <- apply(samples, 1, function(p) {
## 	r_g <- dat$state2 - with(dat, 1 - (1 - predict.stm_point(p[1:5], annual_mean_temp, tot_annual_pp))^(interval/5))
## 	r_e <- as.integer(!(dat$state2)) - 
## 			with(dat, 1 - (1 - predict.stm_point(p[6:10], annual_mean_temp, tot_annual_pp))^(interval/5))
## 	r <- r_g
## 	r[dat$state1 == 1] <- r_e[dat$state1 == 1]
## 	r
## })
## 
## morI <- function(lag0, lag1, coords, resids) {
## 	r.nb <- dnearneigh(coords, lag0, lag1)
## 	r.lw <- nb2listw(r.nb, nbdists(r.nb, coords), zero.policy=TRUE)
## 	mt <- moran.test(resids, r.lw, randomisation=FALSE, zero.policy=TRUE)
## 	mt$estimate
## }
## 
## # get mean moran's i for each lag
## xx <- lags[2:length(lags)]/1000
## residMeans <- rowMeans(resids)
## coords <- as.matrix(dat[,c('x', 'y')])
## mi <- t(mcmapply(morI, lags[1:(length(lags) - 1)], lags[2:length(lags)], 
## 				MoreArgs = list(coords=coords, resids=residMeans), mc.cores=numCores))
## saveRDS(data.frame(mi, lag=xx), paste0("res/moran/", spName, "_moranMean.rds"))
## 
## # get for some posterior samples
## mi.samp <- apply(resids[,sample(ncol(resids), numIters)], 2, function(r) {
## 	t(mcmapply(morI, lags[1:(length(lags) - 1)], lags[2:length(lags)], 
## 				MoreArgs = list(coords=coords, resids=r), mc.cores=numCores)[1,])
## })
## saveRDS(cbind(xx, mi.samp), paste0("res/moran/", spName, "_moranSamp.rds"))
## 
## ymax = min(1, round(max(mi[,1] + 2*sqrt(mi[,3])) + 0.05, 1))
## pdf(w=6, h=6, file=paste0("img/moran/", spName, ".pdf"))
## 	plot(xx, mi[,3], xlab="Lag (km)", ylab="Moran's I", pch=16, ylim=c(-ymax, ymax), xlim=c(0, max(xx)), main=spName, bty='n')
## 	abline(h=0, lty=2)
## 	segments(xx, mi[,1] - 2*sqrt(mi[,3]), xx, mi[,1] + 2*sqrt(mi[,3]))
## 	points(xx, mi[,2], pch=20, cex=0.8, col="#3355ff")
## dev.off()