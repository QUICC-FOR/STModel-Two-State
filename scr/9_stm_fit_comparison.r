#!/usr/bin/Rscript
library(coda)
library(biomod2) # for TSS and ROC

spList = readRDS('dat/speciesList.rds')
spList = c('18032-ABI-BAL', '19290-QUE-ALB', '19489-BET-PAP', 'NA-QUE-PRI')
speciesInfo = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE, colClasses='character')
clim = readRDS(file.path('dat', '/plotClimate_scaled.rds'))

compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}
compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}
point_prob = function(pos, e1, e2, d)
{
	# note that pos.list should be an mcmc list with equal-sized chunks of parallel chains
	# e1 and e2 should be scalars - we are doing this just for a single point on the grid
	
	parBase = rep(0, length(d))
	
	md = as.data.frame(t(sapply(1:nrow(pos), function(i)
	{
		params = parBase
		if(sum(d) != length(pos[i,])) 
		{
			print(paste('design string:', paste(d, collapse='')))
			print(paste('doesnt match posterior:'))
			print(pos)
			cat('\n\n')
		}
		params[d == 1] = pos[i,]
		c(c=compute_c(params, e1, e2),
			e=compute_e(params, e1, e2))
	})))
	
	md$lam = md$c - md$e
	
	sum(md$lam > 0)/nrow(md)
}
postcor = function(params, d, e1, e2, pres) {
	parBase = rep(0, length(d))
	p = parBase
	p[d == 1] = params
	md = data.frame(c = compute_c(p, e1, e2), e=compute_e(p, e1, e2))
	md$lam = md$c - md$e
	md$pres = as.integer(md$lam > 0)
	# if all are the same, the SD will be zero, so we flip one at random to get a value instead of NA
	if(length(table(md$pres)) == 1)
	{
		ind = sample(nrow(md), 1)
		md$pres[ind] = as.integer(!md$pres[ind])
	}
	cor(pres, md$pres)
}	


compute_validation_stats = function(spName)
{
	tryCatch({
		baseDir = file.path('species', spName)
		validationDataset = readRDS(file.path(baseDir, 'dat', paste(spName, 'validationData.rds', sep='_')))
		validationDataset = merge(validationDataset, clim)
		calibrationDataset = readRDS(file.path(baseDir, 'dat', paste(spName, 'presence.rds', sep='_')))
		calibrationDataset = merge(calibrationDataset, clim)
		posterior = readRDS(file.path(baseDir, 'res', paste(spName, 'posterior_thinned.rds', sep='_')))
		pos = do.call(rbind, posterior)
		
		## TEMP - small number of rows and plots
		validationDataset = validationDataset[sample(nrow(validationDataset), 100),]
		calibrationDataset = calibrationDataset[sample(nrow(calibrationDataset), 100),]
		pos = pos[sample(nrow(pos), 100),]
		
		sdm = readRDS(file.path(baseDir, 'res', paste(spName, 'sdm.rds', sep='_')))
		spInfo = speciesInfo[speciesInfo$spName == spName,]
		env1 = spInfo$env1
		env2 = spInfo$env2
		design = spInfo$design
		design = sapply(1:nchar(design), function(i) as.integer(substr(design, i, i)))

		# compute suitability scores for the stm and sdm for the calibration nad validation sets
		# all suitability scores are multiplied by 1000 to prepare them for biomod
		plotProbs = 1000*(sapply(1:nrow(validationDataset), function(i) point_prob(pos, validationDataset[i,env1], validationDataset[i,env2], design)))
		calibPlotProbs = 1000*(sapply(1:nrow(calibrationDataset), function(i) point_prob(pos, calibrationDataset[i,env1], calibrationDataset[i,env2], design)))
		sdmPlotProbs = 1000*(predict(sdm, newdata=validationDataset, type='prob')[,2])
		sdmCalibProbs = 1000*(predict(sdm, newdata=calibrationDataset, type='prob')[,2])

		
		# now measure the correlation between the observed and expected for each posterior replicate
		correlations = sapply(1:nrow(pos), function(i) postcor(pos[i,], design, 
				validationDataset[,env1], validationDataset[,env2], validationDataset[, 3]))

		### see script 4 -- need to do it on calibration data first to determine the threshold
		roc.stm.calib = Find.Optim.Stat(Stat="ROC", Fit=calibPlotProbs, Obs=calibrationDataset[,3])
		roc.stm.valid = Find.Optim.Stat(Stat="ROC", Fit=plotProbs, Obs=validationDataset[,3], Fixed.thresh=roc.stm.calib[2])
		roc.sdm.calib = Find.Optim.Stat(Stat="ROC", Fit=sdmCalibProbs, Obs=calibrationDataset[,3])
		roc.sdm.valid = Find.Optim.Stat(Stat="ROC", Fit=sdmPlotProbs, Obs=validationDataset[,3], Fixed.thresh=roc.sdm.calib[2])


		tss.stm.calib = Find.Optim.Stat(Stat="TSS", Fit=calibPlotProbs, Obs=calibrationDataset[,3])
		tss.stm.valid = Find.Optim.Stat(Stat="TSS", Fit=plotProbs, Obs=validationDataset[,3], Fixed.thresh=tss.stm.calib[2])
		tss.sdm.calib = Find.Optim.Stat(Stat="TSS", Fit=sdmCalibProbs, Obs=calibrationDataset[,3])
		tss.sdm.valid = Find.Optim.Stat(Stat="TSS", Fit=sdmPlotProbs, Obs=validationDataset[,3], Fixed.thresh=tss.sdm.calib[2])

		list(spName=spName, stats=data.frame(stm.roc=roc.stm.valid[1], 
				stm.tss=tss.stm.valid[1], sdm.roc=roc.sdm.valid[1], 
				sdm.tss=tss.sdm.valid[1]), cor=correlations)
				
	}, error=function(e) NA)
}

vStats = lapply(spList, compute_validation_stats)
inds = which(!is.na(vStats))
vStats = vStats[inds]

corStats = sapply(vStats, function(x) x$cor)
colnames(corStats) = spList[inds]

fitStats = as.data.frame(t(sapply(vStats, function(x) x$stats)))
fitStats$species = spList[inds]

dir.create('res')
saveRDS(fitStats, 'res/stmFitComparison.rds')
saveRDS(corStats, 'res/stmFitCorrelations.rds')



	





