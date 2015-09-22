#!/usr/bin/Rscript
library(rgdal)
library(raster)

spInfo = read.csv("dat/speciesInfo.csv", stringsAsFactors=FALSE)
spList = readRDS('dat/speciesList.rds')

compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}


fit_global_models = function(resp, colPredictors, extPredictors, x.coord, y.coord)
{
	require(gam)
	# fit c and e models without the prevalence
	# fits normal (non-spatial) models as well as spatial models
	# resp should be a 2-column matrix (state1 and state2)
	c.sel = which(resp[,1] == 0)
	e.sel = which(resp[,1] == 1)
	
	col.df = data.frame(y=as.integer(resp[c.sel,2] == 1), colPredictors[c.sel,])
	ext.df = data.frame(y=as.integer(resp[e.sel,2] == 0), extPredictors[e.sel,])
	
	c.mod = glm(y ~ ., data=col.df, family=binomial)
	c.mod.xy = gam(y ~ . + s(x.coord[c.sel]) + s(y.coord[c.sel]), 
			data=col.df, family=binomial)
	e.mod = glm(y ~ ., data=ext.df, family=binomial)
	e.mod.xy = gam(y ~ . + s(x.coord[e.sel]) + s(y.coord[e.sel]), 
			data=ext.df, family=binomial)
			
	list(c.mod, c.mod.xy, e.mod, e.mod.xy)
	
}

make_design = function(x, y) paste(paste(x, collapse=""), paste(y, collapse=""), sep="")
parse_design = function(designStr) as.numeric(unlist(strsplit(designStr, "")))

rank_models = function(mods)
{
	result = do.call(rbind, lapply(mods, function(x) data.frame(id = x$id, AIC = x$AIC, 
			BIC = x$BIC, env1 = x$env1, env2 = x$env2, 
			design = make_design(x$colDesign, x$extDesign), stringsAsFactors=FALSE)))
	result = within(result, {
		dAIC = AIC - min(AIC)
		dBIC = BIC - min(BIC)})
	relLL = exp(-0.5 * result$dAIC)
	result$AICw = round(relLL/sum(relLL),4)
	relLL = exp(-0.5 * result$dBIC)
	result$BICw = round(relLL/sum(relLL),4)
	result[order(result$dAIC),]
}

sum_weights = function(ranks)
{
	AICcw = rep(NA, nrow(ranks))
	i = delta = 0
	while(delta <= 10 & i <= nrow(ranks))
	{
		i = i+1
		delta = ranks$dAIC[i]
		if(is.na(AICcw[i]))
		{
			env1 = ranks$env1[i]
			env2 = ranks$env2[i]
			rows = which(ranks$env1 == env1 & ranks$env2 == env2)
			AICcw[rows] = sum(ranks$AICw[rows])
		}
	}
	AICcw
}


select_model = function(ranks, mods, species, useCat=TRUE)
{
	modID = spInfo$modNumber[spInfo$spName == species]
	method = "input from dat/speciesInfo.csv"
	if(modID == 0)
	{
		# find the most complex model where:
		#   1. the cumulative weight is greatest
		#   2. deltaAIC < 2
		method = "automatic selection"
		maxwt = max(ranks$wt_sum, na.rm=TRUE)
		rows = which(ranks$wt_sum == maxwt & ranks$dAIC <= 2)
		complexity = sapply(lapply(ranks$design, parse_design), sum)
		# guard against a model with no predictors, only intercepts
		if(length(rows) == 1 && complexity[rows] == 2)
		{
			maxwt = sort(unique(ranks$wt_sum), TRUE)[2]
			method=paste(method, "; warning: best model had no predictors", sep="")
		}
		rows = which(ranks$wt_sum == maxwt & ranks$dAIC <= 2)		
		if(length(rows) == 0) rows = which(ranks$wt_sum == maxwt)[1]
		modID = ranks$id[rows][which(complexity[rows] == max(complexity[rows]))]
		if(length(modID) > 1) modID = modID[1]
	}
	if(useCat)
	{
		cat("\n")
		print(head(ranks))
		cat(paste("\n selected model", modID, "by", method, "\n"))
	}
	# find the proper model
	mods[[which(lapply(mods, function(x) x$id) == modID)]]
}


prep_species = function(spName)
{
	baseDir = file.path('species', spName)
	annealDir = file.path(baseDir, 'res', 'anneal')
	modelFiles = list.files(annealDir)

	models = lapply(modelFiles, function(x) readRDS(file.path(annealDir, x)))
	models = lapply(1:length(models), function(i) c(models[[i]], list(id=i)))
	

	# set up a data frame with all models ranked
	mRanks = rank_models(models)
	
	# compute the cumulative ranking for all models with the same vars
	mRanks$wt_sum = sum_weights(mRanks)
	
	# select a model
	theMod = select_model(mRanks, models, spName)
		
	saveRDS(mRanks, file.path(baseDir, 'res', paste(spName, 'modelSelection.rds', sep='_')))
	saveRDS(theMod, file.path(baseDir, 'res', paste(spName, 'bestAnnealModel.rds', sep='_')))

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
	params[[spName]] = rep(0, length(design))
	params[[spName]][design == 1] = pars
	env_vars[[spName]]$env1 = theMod$env1
	env_vars[[spName]]$env2 = theMod$env2


	# make global models using no prevalence
	make_preds = function(x1, x2, design)
	{
		vars = data.frame(env1=x1, env2=x2, env12=x1^2, env22=x2^2, env13=x1^3, env23=x2^3)
		design = design[-1]
		vars[, -(which(design == 0))]
	}
	gmods = fit_global_models(resp = cbind(mcmcData$initial, mcmcData$final),
			colPredictors = make_preds(mcmcData$env1, mcmcData$env2, theMod$colDesign),
			extPredictors = make_preds(mcmcData$env1, mcmcData$env2, theMod$extDesign),
			x.coord = transitions$lon, y.coord = transitions$lat)
	saveRDS(gmods, file.path(baseDir, 'res', paste(spName, 'spatialModels.rds', sep='_')))


	mcmcInits = data.frame(
		name = parNames,
		initialValue = params[[spName]],
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = as.integer(!design))
	mcmcInits$priorSD[which(substr(parNames, nchar(parNames), nchar(parNames)) == '0')] = 10

	mcmcDataFile = file.path(baseDir, 'dat', 'mcmc_data.txt')
	mcmcInitFile = file.path(baseDir, 'dat', 'mcmc_inits.txt')

	write.csv(mcmcData, mcmcDataFile, row.names=FALSE)
	write.csv(mcmcInits, mcmcInitFile, row.names=FALSE)
}

env_vars = params = list()
for(species in spList)
{
	cat(paste("Starting species", species, "\n"))
	err_sp = function(e, species) cat(paste("  An error occurred processing species", species, "\n", e, "\n"))
	tryCatch(prep_species(species), error=function(e) err_sp(e, species))
}



# make some prediction maps and response curves
nr = 5
nc = 8

env1.yc = env1.ye = env2.yc = env2.ye = list()
xx = seq(-3, 3, 0.01)
pdf("img/anneal_env1_response.pdf", w=nc*4.5, h=nr*5)
par(mfrow=c(nr, nc), oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))

for(spName in spList)
{
	cat(paste("Plotting species", spName, "\n"))
	p = params[[spName]]
	# find the point at which env1 maximizes lambda
	env1.yc[[spName]] = compute_c(p, xx, rep(0, length(xx)))
	env1.ye[[spName]] = compute_e(p, xx, rep(0, length(xx)))
	env1.lam = env1.yc[[spName]] - env1.ye[[spName]]
	env1.maxx = xx[which(env1.lam == max(env1.lam))[1]]

	# find the point at which env2 maximizes lambda with the max from the prev step
	env2.yc[[spName]] = compute_c(p, rep(env1.maxx, length(xx)), xx)
	env2.ye[[spName]] = compute_e(p, rep(env1.maxx, length(xx)), xx)
	env2.lam = env2.yc[[spName]] - env2.ye[[spName]]
	env2.maxx = xx[which(env2.lam == max(env2.lam))[1]]

	# re-do the computation of xx with the new max from xx
	env1.yc[[spName]] = compute_c(p, xx, rep(env2.maxx, length(xx)))
	env1.ye[[spName]] = compute_e(p, xx, rep(env2.maxx, length(xx)))

	ylim1=range(c(env1.yc[[spName]], env1.ye[[spName]]))
	plot(xx, env1.yc[[spName]], col='blue', xlab=env_vars[[spName]]$env1, ylab="pr", type='l', ylim=ylim1, main=spName)
	lines(xx, env1.ye[[spName]], col='red')
}
dev.off()	

pdf("img/anneal_env2_response.pdf", w=nc*4.5, h=nr*5)
par(mfrow=c(nr, nc), oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))

for(spName in spList)
{
	ylim2=range(c(env2.yc[[spName]], env2.ye[[spName]]))
	plot(xx, env2.yc[[spName]], col='blue', xlab=env_vars[[spName]]$env2, ylab="pr", type='l', ylim=ylim2, main=spName)
	lines(xx, env2.ye[[spName]], col='red')	
}
dev.off()



# map
climDat = readRDS("dat/climateGrid_scaled.rds")
load('dat/map_projections.rdata')
env1 = climDat[, theMod$env1]
env2 = climDat[, theMod$env2]
dpi=600
lamColors = colorRampPalette(c("#008837", "#a6dba0", "#f7f7f7", "#c2a5cf", "#7b3294"), interpolate='spline', space="rgb", bias=1.0)(200)
png(width=as.integer(dpi*4.5*nc), height=as.integer(dpi*5*nr), file="img/anneal_response_map.png", pointsize=12, res=dpi)
par(mfrow=c(nr, nc), oma=c(0,0,2,0), mar=c(4.5,4.5,0.5,0.5))

for(spName in spList)
{
	p = params[[spName]]
	lamVals = compute_c(p, env1, env2) - compute_e(p, env1, env2)
	lambda = data.frame(lon = climDat$lon, lat = climDat$lat, lambda=lamVals)
	at = pretty(range(lamVals))
	labels = as.character(at)
	arg = list(at=at, labels=labels)
	breaks = c(seq(min(lamVals), 0, length.out=100), seq(0.001, max(lamVals), length.out=100))
	
	coordinates(lambda) = c('lon', 'lat')
	gridded(lambda) = TRUE
	lambda = raster(lambda)
	proj4string(lambda) = P4S.latlon
	lambda = projectRaster(lambda, crs=stmMapProjection)
	plot(lambda, col=lamColors, xaxt='n', yaxt='n', main=spName, breaks=breaks, axis.arg=arg, legend=FALSE)
}
dev.off()
