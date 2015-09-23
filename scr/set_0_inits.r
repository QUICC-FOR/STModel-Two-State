#!/usr/bin/Rscript
for(spName in spList)
{
	baseDir = file.path('species', spName)
	theMod = readRDS(file.path(baseDir, 'res', paste(spName, 
			'bestAnnealModel.rds', sep='_')))

	parNames = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
	design =c(theMod$colDesign, theMod$extDesign)
	allParams = rep(0, length(design))

	mcmcInits = data.frame(
		name = parNames,
		initialValue = 0,
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = as.integer(!design))
	mcmcInits$priorSD[which(substr(parNames, nchar(parNames), nchar(parNames)) == '0')] = 10
	mcmcInitFile = file.path(baseDir, 'dat', 'mcmc_inits.txt')
	write.csv(mcmcInits, mcmcInitFile, row.names=FALSE)
}