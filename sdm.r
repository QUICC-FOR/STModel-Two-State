#!/usr/bin/Rscript

inputFile = "dat/statesTwoState.csv"
library(gam)
if(interactive()) {
	# clean workspace to make sure the script is self-contained
	graphics.off()
	rm(list=ls()[ls() != "inputFile"])
	setwd("/Users/mtalluto/Documents/git_projects/Two-State-Metapopulation-Range-Model")
} else {
	arg = commandArgs(TRUE)
	if(length(arg) > 0)
		inputFile = arg[1]
}

dat = read.csv(inputFile, row.names=NULL)
dat.scaled = cbind(dat[,1:3], sapply(dat[,4:7], scale))

# data exploration; only performed for interactive sessions
# if(interactive()) {
# 	par(mfrow=c(2,2))
# 	plot(dat$meanTemp, dat$state, col=(1+dat$state), pch=16, cex=0.7)
# 	plot(dat$maxTemp, dat$state, col=(1+dat$state), pch=16, cex=0.7)
# 	plot(dat$minTemp, dat$state, col=(1+dat$state), pch=16, cex=0.7)
# 	plot(dat$annualPrecip, dat$state, col=(1+dat$state), pch=16, cex=0.7)
# 	cor(dat[,4:7])
# }

sdm.glm = step(glm(state ~ annualPrecip + I(annualPrecip^2) + I(annualPrecip^3) + meanTemp + 
		I(meanTemp^2) + I(meanTemp^3), data=dat.scaled, family=binomial), k=log(nrow(dat.scaled)))
		
sdm.gam = step.gam(gam(state ~ s(annualPrecip) + s(meanTemp), data=dat.scaled, family=binomial),
		scope = list(~1+annualPrecip + s(annualPrecip), ~1+meanTemp+s(meanTemp)))

if(interactive()) {
	# plot response curves
	xx = seq(-3,3, 0.05)
	d1 = data.frame(meanTemp=xx, annualPrecip=rep(median(dat.scaled$annualPrecip[dat.scaled$state==1]), length(xx)))
	par(mfrow=c(1,2), pch=16, cex=0.6)
	plot(dat.scaled$meanTemp, dat.scaled$state, col=(1+dat.scaled$state), xlab="Mean Temp", ylab="pr(presence)")
	lines(d1$meanTemp, predict(sdm.glm, type='response', newdata=d1))
	lines(d1$meanTemp, predict(sdm.gam, type='response', newdata=d1), col='blue')
	
	d1 = data.frame(meanTemp=rep(median(dat.scaled$meanTemp[dat.scaled$state==1]), length(xx)), annualPrecip=xx)
	plot(dat.scaled$annualPrecip, dat.scaled$state, col=(1+dat.scaled$state), xlab="Annual Precip", ylab="pr(presence)")
	lines(d1$annualPrecip, predict(sdm.glm, type='response', newdata=d1))
	lines(d1$annualPrecip, predict(sdm.gam, type='response', newdata=d1),col='blue')
	
}

save(sdm.glm, sdm.gam, file="sdm.rdata")