#!/usr/bin/env Rscript
## library(spaMM)
library(coda)
library(spaMM)
## build spatial autoregressive (and nonspatial) versions of the STM to compare parameter estimates

# needed files and directories
# dat/stm_calib/*

if(interactive()) {
	setwd("/Users/mtalluto/Documents/work/projects_active/stm/STModel-Two-State")
}

# all species > 0.25
spList = c('32945-FRA-NIG', '19287-QUE-MAC', '183397-TSU-CAN', '183375-PIN-RES', '22463-POP-GRA', '183412-LAR-LAR')


# for testing
## spName <- spList[1]

dat_sel <- function(dat, ratio = 2) {
	st1 <- dat[1, 'state1']
	st2 = as.integer(!(st1))
	trans <- which(dat$state2 == st2)
	notrans <- which(dat$state2 == st1)
	nsamp <- as.integer(ratio * length(trans))
	notrans_sel <- sample(notrans, nsamp)
	dat[c(trans, notrans_sel),]
}

mod_plot <- function(mod, mod2, mod3, sep=0.1, ...) {
	coefs <- summary(mod)$coefficients[,1:2]
	coefs2 <- summary(mod2)$coefficients[,1:2]
	coefs3 <- mod3
	allc <- rbind(coefs, coefs2, coefs3)
	yl <- c(min(allc[,1] - allc[,2]), max(allc[,1] + allc[,2]))
	plot(1:5, coefs[,1], pch=20, xlab="Parameter", ylim=yl, col='red', ...)
	segments(1:5, coefs[,1] - coefs[,2], 1:5, coefs[,1] + coefs[,2], col='red')
	points(1:5+sep, coefs2[,1], pch=20, col='blue')
	segments(1:5+sep, coefs2[,1] - coefs2[,2], 1:5+sep, coefs2[,1] + coefs2[,2], col='blue')
	points(1:5+sep+sep, coefs3[,1], pch=20, col='cyan')
	segments(1:5+sep+sep, coefs3[,1] - coefs3[,2], 1:5+sep+sep, coefs3[,1] + coefs3[,2], col='cyan')
	
}

morI <- function(lag0, lag1, coords, resids) {
	r.nb <- dnearneigh(coords, lag0, lag1)
	r.lw <- nb2listw(r.nb, nbdists(r.nb, coords), zero.policy=TRUE)
	mt <- moran.test(resids, r.lw, randomisation=FALSE, zero.policy=TRUE)
	mt$estimate
}

spName = spList[3]


resp_curve <- function(mod1, mod2, ...) {
	xx <- data.frame(annual_mean_temp=seq(-3, 3, length.out=1000), tot_annual_pp=0)
	y1 <- predict(mod1, newdata=xx, type='response')
	y2 <- predict(mod2, newdata=xx, type='response')
	plot(xx[,1], y1, type='l', xlab="Temperature", col='red', ylim=c(0,1), ...)
	lines(xx[,1], y2, col='blue')
}

lags <- 1000 * c(0,10)
par(mfcol=c(4,6))
sapply(spList, function(spName) {
	# first subsample data and compare the subsampled to an full regression
	calib <- readRDS(paste0('dat/stm_calib/', spName, '_stm_calib.rds'))
	calib$interval <- calib$year2 - calib$year1
	samp <- readRDS(paste0('res/posterior/', spName, '_0_samples.rds'))
	colMCMC <- cbind(colMeans(samp[,1:5]), apply(samp[,1:5], 2, sd))
	extMCMC <- cbind(colMeans(samp[,6:10]), apply(samp[,6:10], 2, sd))
	colDat <- calib[calib$state1 == 0,]
	colDat$transition <- colDat$state2 # a 'success' is a transition, which is when state2 is 1
	colDat_sel <- dat_sel(colDat, 10)
	extDat <- calib[calib$state1 == 1,]
	extDat$transition <- as.integer(!(extDat$state2)) # here a success an extinction, which is when state2 is a 0
	extDat_sel <- dat_sel(extDat)

	# fit the 4 models
	form <- "transition ~ prevalence + interval + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + I(tot_annual_pp^2)"
	form_e <- "transition ~ interval + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + I(tot_annual_pp^2)"
	colMod <- glm(as.formula(form), data=colDat, family=binomial)
	colMod_sel <- glm(as.formula(form), data=colDat_sel, family=binomial)
	extMod <- glm(as.formula(form_e), data=extDat, family=binomial)
	extMod_sel <- glm(as.formula(form_e), data=extDat_sel, family=binomial)

## 	mod_plot(colMod, colMod_sel, colMCMC, ylab='Colonization')
## 	mod_plot(extMod, extMod_sel, extMCMC, ylab='Extinction')
## 	resp_curve(colMod, colMod_sel, ylab='Colonization Prob')
## 	resp_curve(extMod, extMod_sel, ylab='Extinction Prob')

	# fit spatial models
	coords_sel <- as.matrix(colDat_sel[,c('x', 'y')])
	r.nb_sel <- dnearneigh(coords_sel, lags[1], lags[2])
	r.lw_sel <- nb2listw(r.nb_sel, nbdists(r.nb_sel, coords_sel), zero.policy=TRUE)
	colMod_sp_sel <- HLCor(transition ~ prevalence + interval + annual_mean_temp + tot_annual_pp + 
					I(annual_mean_temp^2) + I(tot_annual_pp^2) + Matern(1|x+y), data=colDat_sel, family=binomial(), 
					ranPars=list(nu=0.5, rho=1/0.7))


})

### so the answer is that subsamplign changes the intercept. it sucks, but what are you gonna do?
## in the end, I will still do it, because it's the only way to fit the models. the idea is to examine whether the spatial
## pattern matters; maybe it's not such a big deal if the intercept parameter is biased



















## old stuff below
















## source("scr/stm_functions.r")
lags <- 1000 * c(0,200)
datSampleSize <- 2000

# get data
datAll <- readRDS(paste0("dat/stm_calib/", spName, "_stm_calib.rds"))
datAll$interval = datAll$year2 - datAll$year1

# subsample data
dat <- datAll[sample(nrow(datAll), datSampleSize),]

dat_c <- dat[dat$state1 == 0,]

coords <- as.matrix(dat_c[,c('x', 'y')])
r.nb <- dnearneigh(coords, lags[1], lags[2])
r.lw <- nb2listw(r.nb, nbdists(r.nb, coords), zero.policy=TRUE)
lagsarlm(state2 ~ 

spam <- HLCor(state2 ~ annual_mean_temp + tot_annual_pp + interval + prevalence + Matern(1|x+y), data = dat_c, 
	family=binomial(), ranPars=list(nu=0.5, rho=1/0.7))

# nonspatial model for comparison
nospam <- glm(state2 ~ annual_mean_temp + tot_annual_pp + interval + prevalence, data = dat_c, family=binomial)
