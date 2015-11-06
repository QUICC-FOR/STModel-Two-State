#!/usr/bin/env Rscript

# compares posterior response curves to CE rates using AUC & TSS
library(coda)
library(foreach)
source('scr/stm_functions.r')
speciesList = readRDS("dat/speciesList.rds")
## STEPS
#
# foreach species
# √√1. read in calib and valid datasets and posteriors
# foreach model type
#    foreach posterior replicate
#       √√2. compute colonization prob (remembering prevalence, if needed)
#       √√3. compute extinction prob
#       4. multiply suitability scores by 1000
#       5. compute AUC and TSS on calibration data
#       6. compute AUC and TSS on validation data, using the calib results for the threshold
#




interval_transform = function(probs, destination, current) 1 - (1 - probs)^(destination/current)

ceprobs = function(dat, pars, env1=1, env2=2, prev=3, interval=4,
		prevalenceType=c('sdm', 'stm', NA), parInterval=5)
{
	# computes colonization and extinction probs for a dataset given a set of parameters
	# dat: the data, including env1, env2, prevalence, and the time interval
	# env1, env2, prev, interval: column numbers or names in dat
	# prevalenceType: defaults to sdm, which uses the prevalence column in dat
	#   stm will use the equilibrium prevalence, 1 - E/C
	#   NA will use no prevalence
	# parInterval: time interval at which the parameters were estimated
	#
	# value: a data frame with same number of rows as dat and 2 columns, the C and E probs
	
	C.pars = if(length(pars) > 2) pars[1:5] else pars[1]
	E.pars = if(length(pars) > 2) pars[6:10] else pars[2]
	C = predict.stm_point(C.pars, dat[,env1], dat[,env2])
	E = predict.stm_point(E.pars, dat[,env1], dat[,env2])
	Ci = interval_transform(C, dat[,interval], parInterval)
	Ei = interval_transform(E, dat[,interval], parInterval)
	
	if(is.na(prevalenceType)) {
		prevalence = 1 
	} else if(prevalenceType == 'stm') {
		prevalence = 1 - (E/C)
		prevalence[prevalence < 0] = 0
	} else {
		prevalence = dat[,prev]
	}
	C * prevalence
	data.frame(C = Ci * prevalence, E = Ei)
}

auctss = function(calib, valid, initial = 1, final = 2, C = 3, E = 4, biomodScaling = 1000)
{
	# compute auc and tss
	# calib: calibration dataset; must have a column for initial state, final state
	#      and C and E probs
	# valid: validation dataset, same structure (including column names/orders) as calib
	# C, E, initial, final: column names/numbers in the calib and valid dataframes
	
	require(biomod2)
	# we need separate C and E datasets
	calib.c.ind = which(calib[,initial] == 0)
	calib.e.ind = which(calib[,initial] == 1)
## 	valid.c.ind = which(valid[,initial] == 0)
## 	valid.e.ind = which(valid[,initial] == 1)
	
	# compute whether a c/e was observed
	calib$obs = 0
	calib$obs[calib.c.ind][calib[calib.c.ind, final] == 1] = 1
	calib$obs[calib.e.ind][calib[calib.e.ind, final] == 0] = 1
## 	valid$obs = 0
## 	valid$obs[valid.c.ind][valid[valid.c.ind, final] == 1] = 1
## 	valid$obs[valid.e.ind][valid[valid.e.ind, final] == 0] = 1
	
	# get fits
	calib$fit = 0
	calib$fit[calib.c.ind] = calib[calib.c.ind,C]
	calib$fit[calib.e.ind] = calib[calib.e.ind,E]
## 	calib$fit = calib$fit * biomodScaling
## 	valid$fit = 0
## 	valid$fit[valid.c.ind] = valid[valid.c.ind,C]
## 	valid$fit[valid.e.ind] = valid[valid.e.ind,E]
## 	valid$fit = valid$fit * biomodScaling

	print(range(calib$fit))
## 	calib = calib[calib$fit < 200,]
	
	cdat = data.frame(seq(1:nrow(calib)), calib$obs, calib$fit)
	cdat
## 	roc.plot.calculate(cdat)
## 	roc.calib = list(
## 		C = Find.Optim.Stat(Stat="ROC", Fit=calib$fit[calib.c.ind], Obs=calib$obs[calib.c.ind]),
## 		E = Find.Optim.Stat(Stat="ROC", Fit=calib$fit[calib.e.ind], Obs=calib$obs[calib.e.ind]),
## 		all = Find.Optim.Stat(Stat="ROC", Fit=calib$fit, Obs=calib$obs))
## 	roc.valid = list(
## 		C = Find.Optim.Stat(Stat="ROC", Fit=valid$fit[valid.c.ind], Obs=valid$obs[valid.c.ind], Fixed.thresh=roc.calib$C[2]),
## 		E = Find.Optim.Stat(Stat="ROC", Fit=valid$fit[valid.e.ind], Obs=valid$obs[valid.e.ind], Fixed.thresh=roc.calib$E[2]),
## 		all = Find.Optim.Stat(Stat="ROC", Fit=valid$fit, Obs=valid$obs, Fixed.thresh=roc.calib$all[2]))

## 	tss.calib = list(
## 		C = Find.Optim.Stat(Stat="TSS", Fit=calib$fit[calib.c.ind], Obs=calib$obs[calib.c.ind]),
## 		E = Find.Optim.Stat(Stat="TSS", Fit=calib$fit[calib.e.ind], Obs=calib$obs[calib.e.ind]),
## 		all = Find.Optim.Stat(Stat="TSS", Fit=calib$fit, Obs=calib$obs))
## 	tss.valid = list(
## 		C = Find.Optim.Stat(Stat="TSS", Fit=valid$fit[valid.c.ind], Obs=valid$obs[valid.c.ind], Fixed.thresh=roc.calib$C[2]),
## 		E = Find.Optim.Stat(Stat="TSS", Fit=valid$fit[valid.e.ind], Obs=valid$obs[valid.e.ind], Fixed.thresh=roc.calib$E[2]),
## 		all = Find.Optim.Stat(Stat="TSS", Fit=valid$fit, Obs=valid$obs, Fixed.thresh=roc.calib$all[2]))

## 	c(roc.C = roc.calib$C[1], roc.E = roc.calib$E[1], roc.all = roc.calib$all[1], 
## 			tss.C = tss.calib$C[1], tss.E = tss.calib$E[1], tss.all = tss.calib$all[1])
## 	c(roc.C = roc.valid$C[1], roc.E = roc.valid$E[1], roc.all = roc.valid$all[1], 
## 			tss.C = tss.valid$C[1], tss.E = tss.valid$E[1], tss.all = tss.valid$all[1])
}




compute_evaluation = function(calib, valid, params, modName)
{
	# 2,3. compute colonization and extinction probs
	prType = 'sdm'
	if(modName == 'a' | modName == 'ia') prType = 'stm'
	if(modName == 'g' | modName == 'ig') prType = NA
	calib$interval = calib$year2 - calib$year1
	valid$interval = valid$year2 - valid$year1
	calib.CE = ceprobs(calib, params, env1 = 'annual_mean_temp', env2 = 'tot_annual_pp',
			prev = 'prevalence', interval = 'interval',, prevalenceType = prType)
	valid.CE = ceprobs(valid, params, env1 = 'annual_mean_temp', env2 = 'tot_annual_pp',
			prev = 'prevalence', interval = 'interval',, prevalenceType = prType)

	# 4-6. get auc/tss
	cdat = cbind(calib$state1, calib$state2, calib.CE)
	vdat = cbind(valid$state1, valid$state2, valid.CE)
	
	
	auctss(cdat, vdat)
}

spName = speciesList[1]
# 1. read in calib and valid datasets and posterior samples
calib = readRDS(file.path('dat','stm_calib', paste0(spName,'stm_calib.rds'))) 
valid = readRDS(file.path('dat','stm_valid', paste0(spName,'stm_valid.rds'))) 
posterior = readRDS(file.path('res', 'posterior', paste0(spName, '_posterior.rds')))

modName = '0'
samples = posterior[[modName]]
## subsample = 4
## samples = samples[sample(nrow(samples), subsample),]

## calibration curve
cut.mid <- function(x, breaks) {
  r <- range(x)
  b <- seq(r[1], r[2], length=breaks+1)
  bottom = b[1:breaks]
  top = b[2:(breaks+1)]
  mid = rowMeans(cbind(b[1:breaks], b[2:(breaks+1)]))



  k <- cut(x, breaks=breaks)
  levels(k) = mid
  as.numeric(as.character(k))
}


i = 1000
params = samples[i,]
prType = 'sdm'
if(modName == 'a' | modName == 'ia') prType = 'stm'
if(modName == 'g' | modName == 'ig') prType = NA
calib$interval = calib$year2 - calib$year1
calib.CE = ceprobs(calib, params, env1 = 'annual_mean_temp', env2 = 'tot_annual_pp',
		prev = 'prevalence', interval = 'interval',, prevalenceType = prType)
calib.c.ind = which(calib$state1 == 0)
calib.e.ind = which(calib$state1 == 1)
calib$obs = 0
calib$obs[calib.c.ind][calib$state2[calib.c.ind] == 1] = 1
calib$obs[calib.e.ind][calib$state2[calib.e.ind] == 0] = 1

predicted = cut.mid(calib.CE[calib.e.ind,2], 100)
obs = aggregate(calib$obs[calib.e.ind], by=list(pred=predicted), FUN = function(x) sum(x)/length(x))

n = table(predicted)
cex = log(n)/log(max(n)) + .2
cex2 =.6
par(mfrow=c(1,3))
plot(obs$pred, obs$x, pch=16, xlim=c(0, max(obs$pred)), ylim=c(0, max(obs$x)), xlab='predicted', ylab='fit', cex=cex)
abline(a=0, b=1, lty=2)
plot(obs$pred, obs$x, pch=16, xlim=c(0, max(obs$pred)), ylim=c(0, max(obs$x)), xlab='predicted', ylab='fit', cex=cex2)
abline(a=0, b=1, lty=2)

cols = rep('#4444ffcc', length(calib.e.ind))
cols[calib$state2[calib.e.ind] == 0] = '#ff4444ff'
cex3 = rep(0.25, length(calib.e.ind))
cex3[calib$state2[calib.e.ind] == 0] = 0.5
plot(calib$annual_mean_temp[calib.e.ind], calib.CE[calib.e.ind,2], pch=16, cex=cex3, col=cols, xlab='temperature', ylab='extinction')



	cdat = compute_evaluation(calib, valid, samples[i,], modName)
	ac = auc(cdat)
	test = roc.plot.calculate(cdat[cdat[,3] < 0.2,], threshold=1000 )
	
	par(mfrow=c(1,2))
plot(obs$pred, obs$x, pch=16, xlim=c(0, max(obs$pred)), ylim=c(0, max(obs$x)), xlab='predicted', ylab='fit', cex=cex)
abline(a=0, b=1, lty=2)
	plot(1-test[,4],test[,3])
abline(a=0, b=1, lty=2)


## library(doParallel)
## registerDoParallel(cores=4)
## foreach(i = 1:nrow(samples), .combine=rbind, .packages='biomod2') %dopar% {
## 		compute_evaluation(calib, valid, samples[i,], modName)
## }
## 

