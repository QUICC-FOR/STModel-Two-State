#!/usr/bin/env Rscript
library(coda)
library(spdep)
library(sp)
library(parallel)

## compute Moran's I for a given input species at a fixed set of lag distances
# pretty slow due to the large number of pairwise computations required

# needed files and directories
# res/moran
# img/moran
# scr/stm_functions.r
# dat/stm_calib/*
# res/posterior/*

if(interactive()) {
	setwd("/Users/mtalluto/Documents/work/projects_active/stm/STModel-Two-State")
	spName = '18032-ABI-BAL'
} else {
	arg = commandArgs(TRUE)
	spName = arg[1]
}

source("scr/stm_functions.r")
lags <- 1000 * seq(0, 700, 25)
numCores = 16
parSampleSize <- 1000
datSampleSize <- 20000
numIters <- 100

# get data
datAll <- readRDS(paste0("dat/stm_calib/", spName, "_stm_calib.rds"))
datAll$interval = datAll$year2 - datAll$year1

# get par samples and subsample
samples <- readRDS(paste0("res/posterior/", spName, "_0_samples.rds"))
samples <- samples[sample(nrow(samples), parSampleSize),]

# subsample data
dat <- datAll[sample(nrow(datAll), datSampleSize),]

resids <- apply(samples, 1, function(p) {
	r_g <- dat$state2 - with(dat, 1 - (1 - predict.stm_point(p[1:5], annual_mean_temp, tot_annual_pp))^(interval/5))
	r_e <- as.integer(!(dat$state2)) - 
			with(dat, 1 - (1 - predict.stm_point(p[6:10], annual_mean_temp, tot_annual_pp))^(interval/5))
	r <- r_g
	r[dat$state1 == 1] <- r_e[dat$state1 == 1]
	r
})

morI <- function(lag0, lag1, coords, resids) {
	r.nb <- dnearneigh(coords, lag0, lag1)
	r.lw <- nb2listw(r.nb, nbdists(r.nb, coords), zero.policy=TRUE)
	mt <- moran.test(resids, r.lw, randomisation=FALSE, zero.policy=TRUE)
	mt$estimate
}

# get mean moran's i for each lag
xx <- lags[2:length(lags)]/1000
residMeans <- rowMeans(resids)
coords <- as.matrix(dat[,c('x', 'y')])
mi <- t(mcmapply(morI, lags[1:(length(lags) - 1)], lags[2:length(lags)], 
				MoreArgs = list(coords=coords, resids=residMeans), mc.cores=numCores))
saveRDS(data.frame(mi, lag=xx), paste0("res/moran/", spName, "_moranMean.rds"))

# get for some posterior samples
mi.samp <- apply(resids[,sample(ncol(resids), numIters)], 2, function(r) {
	t(mcmapply(morI, lags[1:(length(lags) - 1)], lags[2:length(lags)], 
				MoreArgs = list(coords=coords, resids=r), mc.cores=numCores)[1,])
})
saveRDS(cbind(xx, mi.samp), paste0("res/moran/", spName, "_moranSamp.rds"))

ymax = min(1, round(max(mi[,1] + 2*sqrt(mi[,3])) + 0.05, 1))
pdf(w=6, h=6, file=paste0("img/moran/", spName, ".pdf"))
	plot(xx, mi[,3], xlab="Lag (km)", ylab="Moran's I", pch=16, ylim=c(-ymax, ymax), xlim=c(0, max(xx)), main=spName, bty='n')
	abline(h=0, lty=2)
	segments(xx, mi[,1] - 2*sqrt(mi[,3]), xx, mi[,1] + 2*sqrt(mi[,3]))
	points(xx, mi[,2], pch=20, cex=0.8, col="#3355ff")
dev.off()