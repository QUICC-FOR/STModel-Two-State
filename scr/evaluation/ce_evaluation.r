#!/usr/bin/env Rscript

# compares posterior response curves to CE rates using AUC & TSS
library(coda)
library(foreach)
library(doParallel)
source('scr/stm_functions.r')

clArgs <- commandArgs(trailingOnly = TRUE)

## these will be read in from the command line
spName = clArgs[1]
modName = clArgs[2]
numCores = if(length(clArgs) > 2) clArgs[3] else detectCores()
sampleSize = if(length(clArgs) > 3) clArgs[4] else NA

registerDoParallel(cores=numCores)


interval_transform = function(probs, destination, current) 1 - (1 - probs)^(destination/current)

ceprobs = function(dat, pars, env1=1, env2=2, prev=3, interval=4,
		prevalenceType=c('sdm', 'stm', NA), parInterval=5, E = NA)
{
	# computes colonization or extinction probs for a dataset given a set of parameters
	# dat: the data, including env1, env2, prevalence, and the time interval
	# env1, env2, prev, interval: column numbers or names in dat
	# prevalenceType: defaults to sdm, which uses the prevalence column in dat
	#   stm will use the equilibrium prevalence, 1 - E/C
	#   NA will use no prevalence
	# parInterval: time interval at which the parameters were estimated
	#
	# value: a data frame with same number of rows as dat and 2 columns, the C and E probs
	
	vals = predict.stm_point(pars, dat[,env1], dat[,env2])
	Vi = interval_transform(vals, dat[,interval], parInterval)
	
	if(is.na(prevalenceType)) {
		prevalence = 1 
	} else if(prevalenceType == 'stm') {
		prevalence = 1 - (E/vals)
		prevalence[prevalence < 0] = 0
	} else {
		prevalence = dat[,prev]
	}
	vals * prevalence
}


auctss = function(dat, threshold = NULL, biomodScaling = 1000)
{
	# compute auc and tss
	require(biomod2)
	dat$fit = dat$fit * biomodScaling
	roc = Find.Optim.Stat(Stat="ROC", Fit=dat$fit, Obs=dat$obs, Fixed.thresh=threshold[1])
	tss = Find.Optim.Stat(Stat="TSS", Fit=dat$fit, Obs=dat$obs, Fixed.thresh=threshold[2])
	
	c(roc = roc[1], tss = tss[1], roc.thresh = roc[2], tss.thresh = tss[2])
}



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


prep_dat = function(dat, params, modName)
{
	prType = 'sdm'
	if(modName == 'a' | modName == 'ia') prType = 'stm'
	if(modName == 'g' | modName == 'ig') prType = NA

	# compute sampling intervals
	dat = lapply(dat, function(x) within(x, interval <- year2 - year1))

	# get c and e probs
	C.pars = if(length(params) > 2) params[1:5] else params[1]
	E.pars = if(length(params) > 2) params[6:10] else params[2]
	dat = lapply(dat, function(x) within(x, {
		E = ceprobs(x, E.pars, env1 = 'annual_mean_temp', env2 = 'tot_annual_pp',
			prev = 'prevalence', interval = 'interval', prevalenceType = NA)
		C = ceprobs(x, C.pars, env1 = 'annual_mean_temp', env2 = 'tot_annual_pp',
## 			prev = 'prevalence', interval = 'interval', prevalenceType = NA, E=E)
 			prev = 'prevalence', interval = 'interval', prevalenceType = prType, E=E)
	}))
	
	# split into separate C and E datasets and compute observed and fit
	dat = lapply(dat, function(x) {
		# find the indices for e and c datasets
		c.ind = which(x$state1 == 0)
		e.ind = which(x$state1 == 1)
		x$obs = x$fit = 0
		x$obs[c.ind][x$state2[c.ind] == 1] = 1
		x$obs[e.ind][x$state2[e.ind] == 0] = 1
	
		# get fits
		x$fit[c.ind] = x$C[c.ind]
		x$fit[e.ind] = x$E[e.ind]

		list(all = x, cdat = x[c.ind,], edat = x[e.ind,])
	})
	
	dat
}

compute_evaluation = function(dat, params, modName, cc.pred.x)
{
	dat = prep_dat(dat, params, modName)
	
	# get auc and tss
	stats = list()
	stats$calib = lapply(dat$calib, function(x) {
		res = compute_cc(x, cc.pred.x)
		res$auctss = auctss(x)
		res
	})
	stats$valid = lapply(names(dat$valid), function(nm) {
		x = dat$valid[[nm]]
		res = compute_cc(x, cc.pred.x)
		res$auctss = auctss(x, stats$calib[[nm]]$auctss[3:4])
		names(res$auctss) = names(stats$calib[[nm]]$auctss)
		res
	})
	
	stats
}


compute_cc = function(x, pred.x=x$pred, ccbins=50)
{
	# compute the calibration curve
	predicted = cut.mid(x$fit, ccbins)
	obs = aggregate(x$obs, by=list(pred=predicted), FUN = function(y) sum(y)/length(y))
	# sort, just to be sure
	obs = obs[order(obs$pred),]
	obs$n = table(predicted)
	cc = lm(x ~ pred, weights = n, data=obs)
	yy = predict(cc, newdata = list(pred = pred.x))
	res = list(cc_stat = c(coefficients(cc), summary(cc)$r.squared), cc_pred = yy)
	names(res$cc_stat) = c('a', 'b', 'r2')
	res
}






# this function recursively combines N lists with similar structures using rbind
# thus, N lists of lists of vectors will return a list of lists of matrices with N rows each
eval_combine = function(...) {
	x1 = list(...)[[1]]
	if(is.list(x1[[1]])) {
		Map(eval_combine, ...)
	} else
		Map(rbind,...)
}



evalData = list(calib = readRDS(file.path('dat','stm_calib', paste0(spName,'stm_calib.rds'))), 
		valid = readRDS(file.path('dat','stm_valid', paste0(spName,'stm_valid.rds'))))
posterior = readRDS(file.path('res', 'posterior', paste0(spName, '_posterior.rds')))[[modName]]

if(is.na(sampleSize)) {
	samples = posterior
} else
	samples = posterior[sample(nrow(posterior), sampleSize),]


# x locations at which to compute the calibration curve
cc.x = seq(0,1,length.out=500)
results = foreach(i = 1:nrow(samples), .combine=eval_combine, .packages='biomod2', .multicombine=TRUE) %dopar% {
		compute_evaluation(evalData, params=samples[i,], modName=modName, cc.pred.x = cc.x)
}
results$cc.x = cc.x
dname = file.path('res', 'ce_eval')
suppressWarnings(dir.create(dname, recursive = TRUE))
fpath = paste0(dname, '/', spName, '_', modName, '_ceEval.rds')
saveRDS(results, fpath)


