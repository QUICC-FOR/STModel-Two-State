#!/usr/bin/Rscript
library(rgdal)
library(raster)
library(biomod2) # for TSS and ROC

compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}


fit_global_models = function(resp, predictors, x.coord, y.coord)
{
	# fit c and e models without the prevalence
	# fits normal (non-spatial) models as well as spatial models
	# resp should be a 2-column matrix (state1 and state2)
	c.sel = which(resp[,1] == 0)
	e.sel = which(resp[,1] == 1)
	
	col.df = data.frame(y=as.integer(resp[c.sel,2] == 1), predictors[c.sel,])
	ext.df = data.frame(y=as.integer(resp[e.sel],2 == 0), predictors[e.sel,])
	
	c.mod = glm(y ~ ., data=col.df, family=binomial)
	c.mod.xy = glm(y ~ . + poly(x.coord[c.sel],3) + poly(y.coord[c.sel], 3), 
			data=col.df, family=binomial)
	e.mod = glm(y ~ ., data=ext.df, family=binomial)
	e.mod.xy = glm(y ~ . + poly(x.coord[e.sel],3) + poly(y.coord[e.sel], 3), 
			data=ext.df, family=binomial)
			
	list(c.mod, c.mod.xy, e.mod, e.mod.xy)
	
}


#spList = readRDS('dat/speciesList.rds')
spList = c('19254-JUG-NIG', '19050-ULM-RUB', '23690-OXY-ARB', '19511-OST-VIR', '32945-FRA-NIG')
for(spName in spList)
{
	cat(paste("Starting species", spName, "\n"))
	baseDir = file.path('species', spName)
	annealDir = file.path(baseDir, 'res', 'anneal')
	modelFiles = list.files(annealDir)

	models = lapply(modelFiles, function(x) readRDS(file.path(annealDir, x)))
	models = lapply(1:length(models), function(i) c(models[[i]], list(id=i)))

	# compute AUC
	calibrationSet = readRDS(file.path(baseDir, 'dat', paste(spName, 'stm_calib.rds', sep='_')))
	validationSet = readRDS(file.path(baseDir, 'dat', paste(spName, 'stm_valid.rds', sep='_')))

	prep_prediction = function(x, p, e1, e2)
	{
		if(e1 == "NA" | is.na(e1)) {
			env1 = rep(0, nrow(x))
		} else {
			env1 = x[,e1]
		}
		if(e2 == "NA" | is.na(e2)) {
			env2 = rep(0, nrow(x))
		} else {
			env2 = x[,e2]
		}
		
		interval = x$year2 - x$year1
		c.rate = compute_c(p, env1, env2)
		e.rate = compute_e(p, env1, env2)
		
		# get the interval rates (duh!)
		# the parameters were estimated for a 5 year interval
		c.rate = 1 - (1 - c.rate)^(interval/5)
		e.rate = 1 - (1 - e.rate)^(interval/5)

		response = x$state2
		prediction = rep(0, length(response))
		
		prediction[x$state1 == 0] = (c.rate * x$prevalence)[x$state1 == 0]
		prediction[x$state1 == 1] = (1 - e.rate)[x$state1 == 1]
		prediction = 1000*prediction # biomod expects data from 1 to 1000
		
		data.frame(response, prediction)
	}

	compute_valid = function(design, params, env1, env2)
	{
		p = rep(0, length(design))
		p[design == 1] = params
		
		calib = prep_prediction(calibrationSet, p, env1, env2)
		valid = prep_prediction(validationSet, p, env1, env2)
	
		roc.calib = Find.Optim.Stat(Stat="ROC", Fit=calib$prediction, Obs=calib$response)
		roc.valid = Find.Optim.Stat(Stat="ROC", Fit=valid$prediction, Obs=valid$response, Fixed.thresh=roc.calib[2])
		tss.calib = Find.Optim.Stat(Stat="TSS", Fit=calib$prediction, Obs=calib$response)
		tss.valid = Find.Optim.Stat(Stat="TSS", Fit=valid$prediction, Obs=valid$response, Fixed.thresh=tss.calib[2])
		
		data.frame(calibROC = roc.calib[1], validROC = roc.valid[1], 
				calibTSS = tss.calib[1], validTSS = tss.valid[1])
	}
	
	make_design = function(x, y) paste(paste(x, collapse=""), paste(y, collapse=""), sep="")
	mRanks = do.call(rbind, lapply(models, function(x) data.frame(id = x$id, AIC = x$AIC, 
			BIC = x$BIC, env1 = x$env1, env2 = x$env2, 
			design = make_design(x$colDesign, x$extDesign), stringsAsFactors=FALSE)))
			
	mRanks = cbind(mRanks, do.call(rbind, lapply(models, function(x)
		compute_valid(c(x$colDesign, x$extDesign), x$annealParams$par, x$env1, x$env2))))
	
	mRanks = within(mRanks, {
		dAIC = AIC - min(AIC)
		dBIC = BIC - min(BIC)})
	relLL = exp(-0.5 * mRanks$dAIC)
	mRanks$AICw = round(relLL/sum(relLL),4)
	relLL = exp(-0.5 * mRanks$dBIC)
	mRanks$BICw = round(relLL/sum(relLL),4)

	mRanks = mRanks[order(mRanks$dAIC),]
	saveRDS(mRanks, file.path(baseDir, 'res', paste(spName, 'modelSelection.rds', sep='_')))
	theMod = models[[mRanks$id[1]]]
	saveRDS(theMod, file.path(baseDir, 'res', paste(spName, 'bestAnnealModel.rds', sep='_')))

	### HERE - need to global models for the best model
	
	
	

	# get mcmc files ready
	transitions = rbind(readRDS(file.path(baseDir, 'dat', paste(spName, 'stm_calib.rds', sep='_'))),
		readRDS(file.path(baseDir, 'dat', paste(spName, 'stm_valid.rds', sep='_'))))

	parNames = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
	mcmcData= data.frame(
		initial = transitions$state1,
		final = transitions$state2,
		env1 = transitions[,theMod$env1],
		env2 = transitions[,theMod$env2],
		interval = transitions$year2-transitions$year1,
		prevalence1 = transitions$prevalence)



	pars = theMod$annealParams$par
	design =c(theMod$colDesign, theMod$extDesign)
	p = rep(0, length(design))
	p[design == 1] = pars


	mcmcInits = data.frame(
		name = parNames,
		initialValue = p,
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = as.integer(!design))
	mcmcInits$priorSD[which(substr(parNames, nchar(parNames), nchar(parNames)) == '0')] = 10

	# second one initialize to 0
	mcmcInits2 = mcmcInits
	mcmcInits2$initialValue = rep(0, nrow(mcmcInits2))

	# third one initialize to mean values for c and e
	mcmcInits3 = mcmcInits2
	mcmcInits3$initialValue[parNames == 'e0'] = 
			qlogis(with(transitions, sum(state1 == 1 & state2 == 0) / sum(state1 == 1)))
	mcmcInits3$initialValue[parNames == 'g0'] = 
			qlogis(with(transitions, sum(state1 == 0 & state2 == 1) / sum(state1 == 0)))

	mcmcDataFile = file.path(baseDir, 'dat', 'mcmc_data.txt')
	mcmcInitFile1 = file.path(baseDir, 'dat', 'mcmc_inits1.txt')
	mcmcInitFile2 = file.path(baseDir, 'dat', 'mcmc_inits2.txt')
	mcmcInitFile3 = file.path(baseDir, 'dat', 'mcmc_inits3.txt')

	write.csv(mcmcData, mcmcDataFile, row.names=FALSE)
	write.csv(mcmcInits, mcmcInitFile1, row.names=FALSE)
	write.csv(mcmcInits2, mcmcInitFile2, row.names=FALSE)
	write.csv(mcmcInits3, mcmcInitFile3, row.names=FALSE)


	## make a prediction map and response curve
	# response curves
	x1 = seq(-3, 3, 0.01)
	x2 = seq(-3, 3, 0.01)

	# find the point at which env1 maximizes lambda
	env1.yc = compute_c(p, x1, rep(0, length(x1)))
	env1.ye = compute_e(p, x1, rep(0, length(x1)))
	env1.lam = env1.yc - env1.ye
	env1.maxx = x1[which(env1.lam == max(env1.lam))[1]]

	# find the point at which env2 maximizes lambda with the max from the prev step
	env2.yc = compute_c(p, rep(env1.maxx, length(x2)), x2)
	env2.ye = compute_e(p, rep(env1.maxx, length(x2)), x2)
	env2.lam = env2.yc - env2.ye
	env2.maxx = x2[which(env2.lam == max(env2.lam))[1]]

	# re-do the computation of x1 with the new max from x2
	env1.yc = compute_c(p, x1, rep(env2.maxx, length(x1)))
	env1.ye = compute_e(p, x1, rep(env2.maxx, length(x1)))


	pdf(w=7,h=4, file=file.path(baseDir, 'img', 'anneal_response.pdf'))
	par(mfrow=c(1,2), oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))
	ylim1=range(c(env1.yc, env1.ye))
	ylim2=range(c(env2.yc, env2.ye))
	plot(x1, env1.yc, col='blue', xlab=theMod$env1, ylab="pr", type='l', ylim=ylim1)
	lines(x1, env1.ye, col='red')
	plot(x2, env2.yc, col='blue', xlab=theMod$env2, ylab="pr", type='l', ylim=ylim2)
	lines(x2, env2.ye, col='red')	
	mtext(spName, outer=T, line=0.5)
	dev.off()


	# map
	climDat = readRDS("dat/climateGrid_scaled.rds")
	load('dat/map_projections.rdata')
	env1 = climDat[, theMod$env1]
	env2 = climDat[, theMod$env2]
	lamVals = compute_c(p, env1, env2) - compute_e(p, env1, env2)
	lambda = data.frame(lon = climDat$lon, lat = climDat$lat, lambda=lamVals)
	coordinates(lambda) = c('lon', 'lat')
	gridded(lambda) = TRUE
	lambda = raster(lambda)
	proj4string(lambda) = P4S.latlon
	lambda = projectRaster(lambda, crs=stmMapProjection)

	lamColors = colorRampPalette(c("#008837", "#a6dba0", "#f7f7f7", "#c2a5cf", "#7b3294"), interpolate='spline', space="rgb", bias=1.0)(200)
	at = pretty(range(lamVals))
	labels = as.character(at)
	arg = list(at=at, labels=labels)
	breaks = c(seq(min(lamVals), 0, length.out=100), seq(0.001, max(lamVals), length.out=100))
	pdf(w=7,h=7, file=file.path(baseDir, 'img', 'anneal_map.pdf'))
	plot(lambda, col=lamColors, xaxt='n', yaxt='n', main=spName, breaks=breaks, axis.arg=arg)
	dev.off()
}
