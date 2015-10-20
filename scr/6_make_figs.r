library(coda)
#speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

speciesList = c('18032-ABI-BAL','28731-ACE-SAC')
climGrid = readRDS(file.path('dat', 'climateGrid_scaled.rds'))
env1 = climGrid$annual_mean_temp
env2 = climGrid$tot_annual_pp
source('scr/stm_functions.r')

predict.stm_point = function(p, env1, env2)
{
	plogis(p[1] + p[2]*env1 + p[3]*env2 + p[4]*env1^2 + p[5]*env2^2)
}

## predict.stm_point_int = function(p, env1, env2)
## {
## 	c(plogis(p[1]), plogis(p[2]))
## }
## 


predict.stm = function(posterior, env1, env2)
{
	predict.func = if(ncol(posterior) == 2)
	{ 
		predict.stm_point_int
	} else {
		predict.stm_point
	}
	preds = array(NA, dim=c(3,nrow(posterior), length(env1)))
	for(i in 1:nrow(posterior))
	{
		preds[1:2,i,] = t(predict.func(posterior[i,], env1, env2))
		preds[3,i,] = preds[i,1] - preds[i,2]		
	}
	res = data.frame(
		c = colMeans(preds[1,,]),
		e = colMeans(preds[2,,]),
		lam = colMeans(preds[3,,]),
		stm = apply(preds[3,,],2, function(x) sum(x > 0))
	)
}


for(spName in speciesList)
{
	# get the posterior
	spPost = readRDS(file.path('res','posterior',paste0(spName,'_posterior.rds')))[['0']]
	# subsample posterior (for now)
	spPost = spPost[sample(nrow(spPost), 1000),]
	spGrid = readRDS(file.path('res','sdm',paste0(spName,'_sdm_projection.rds')))

	# predict col and ext rates for the climate data & for response curves
	pG = pE = pLam = matrix(NA, nrow=nrow(spPost), ncol=length(env1))
	for(i in 1:nrow(spPost))
	{
		pG[i,] = predict.stm_point(spPost[i,1:5], env1, env2)
		pE[i,] = predict.stm_point(spPost[i,6:10], env1, env2)
		pLam[i,] = pG[i,] - pE[i,]
		cat(i,'\n')
	}
	spGrid$stm = apply(pLam, 2, function(x) sum(x > 0, na.rm=TRUE)/length(x))



	## plots here; dragons!


	paperwidth = 10
	dpi = 600
	hToWRatio = 0.4
	width = as.integer(dpi*paperwidth)
	height = as.integer(width * hToWRatio)
	fontsize = 15
	stPres = make_raster(spGrid$stm, spGrid[,1:2], P4S.latlon, stmMapProjection)
	png(w=width, h=height, file=paste0("img/", spName, "_posterior_maps.png"), pointsize=fontsize, res = dpi)
	par(mfrow=c(1,3))
	pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
	plot(stPres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1.00000001))
	plotbg()
	
	sdPres = make_raster(spGrid$sdm, spGrid[,1:2], P4S.latlon, stmMapProjection)
	plot(sdPres, col=pres.colors(100), xaxt='n', yaxt='n', zlim=c(0,1), legend=FALSE)
	plotbg()
	
	sdCutoff = 0.1
	stCutoff = 0.1
	spGrid$sdmPres = as.integer(spGrid$sdm >= sdCutoff)
	spGrid$stmPres = as.integer(spGrid$stm >= stCutoff)
	spGrid$diseq = NA
	spGrid$diseq[spGrid$sdmPres == 1 & spGrid$stmPres == 0] = 2 # contract
	spGrid$diseq[spGrid$sdmPres == 1 & spGrid$stmPres == 1] = 0 # present
	spGrid$diseq[spGrid$sdmPres == 0 & spGrid$stmPres == 1] = 1 # expand
	
	cat.colors = c('#1f78b4', '#b2df8a', '#fb9a99')
	diseq = make_raster(spGrid$diseq, spGrid[,1:2], P4S.latlon, stmMapProjection)
	plot(diseq, col=cat.colors, xaxt='n', yaxt='n', legend=FALSE)
	plotbg()

dev.off()
		
		
	
}

