library(coda)
library(sp)
library(raster)
library(rgdal)

# temporary
spName = '18032-ABI-BAL'
spText = expression(italic(Abies~balsamifera))
env1 = 'annual_mean_temp'
env2 = 'tot_annual_pp'

load('dat/map_projections.rdata')

prj.ras = function(x)
{
	coordinates(x) = 1:2
	gridded(x) = TRUE
	x = raster(x)
	proj4string(x) = P4S.latlon
	projectRaster(x, crs=stmMapProjection)
}

ocean = readOGR(dsn="dat/ne_50m_ocean", layer="ne_50m_ocean")
ocean = spTransform(ocean, stmMapProjection)
lakes = readOGR(dsn="dat/ne_50m_lakes", layer="ne_50m_lakes")
# grab specific lakes
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]
grLakes = spTransform(grLakes, stmMapProjection)

plotbg = function(txt="", rangeMap=T)
{
	plot(ocean, col="white", add=T)
	plot(grLakes, col="white", add=T)
	mtext(txt, adj=0, cex = 1)
	
	# plot grid lines and labels
	llgridlines(ocean, easts=seq(-90, -40, 10), norths=seq(20,50,10), ndiscr=100, 
			side="ES", offset=1e6, col='#44444444')
	par(xpd=TRUE)
#	westText = c("30°N", "40°N", "50°N")
	eastText = c("20°N", "30°N", "40°N")
#	westCoords = list(rep(-3e5,3), c(0.75e6, 1.9e6, 3e6))
	eastCoords = list(rep(3.5e6,3), c(0.55e6, 1.7e6, 2.9e6))
	southText = c("90°W", "80°W", "70°W")
	southCoords = list(c(0.65e6, 1.7e6, 2.8e6), rep(-0.7e6,3))
	cex.ll = 0.8
	text(eastCoords[[1]], eastCoords[[2]], eastText, cex=cex.ll)
#	text(westCoords[[1]], westCoords[[2]], westText, cex=cex.ll)
	text(southCoords[[1]], southCoords[[2]], southText, cex=cex.ll)
	par(xpd=FALSE)
}


compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}


rawPosterior = readRDS(file.path('species', spName, 'res', paste(spName, 'posterior.rds', sep='_')))
posterior = rawPosterior[sample(nrow(rawPosterior), 500),]
climGrid = readRDS('dat/climateGrid_scaled.rds')

posterior.c = sapply(1:nrow(posterior), function(x)
{
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_c(params, climGrid[,env1], climGrid[,env2])
})
posterior.e = sapply(1:nrow(posterior), function(x)
{
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_e(params, climGrid[,env1], climGrid[,env2])
})
posterior.lam = posterior.c - posterior.e

grid.c = data.frame(lon=climGrid$lon, lat=climGrid$lat, c=rowMeans(posterior.c))
grid.e = data.frame(lon=climGrid$lon, lat=climGrid$lat, e=rowMeans(posterior.e))
grid.lam = data.frame(lon=climGrid$lon, lat=climGrid$lat, lam=apply(posterior.lam, 1, function(x) sum(x > 0)/ncol(posterior.lam)))
grid.sdm = readRDS(file.path('species', spName, 'res', paste(spName, 'sdm_grid_projection.rds', sep='_')))
grid.c = prj.ras(grid.c)
grid.e = prj.ras(grid.e)
grid.lam = prj.ras(grid.lam)
grid.sdm = prj.ras(grid.sdm)

pdf(w=6.5, h=8.2, file="img/posterior_maps.pdf")
par(mfrow=c(5,4), mar=c(2,2,2,1), oma=c(0,0,0,2))
cex.title = 0.5
leg.args = list(cex.axis=0.6)
for(i in 1:5)
{
	col.colors = colorRampPalette(c('#eff3ff', '#6baed6', '#08519c'), space="rgb")
	ext.colors = colorRampPalette(c('#fee5d9', '#fb6a4a', '#a50f15'))
	pres.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
			interpolate='spline', bias=2, space="rgb")
	sdm.colors = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#cc99ff"), 
			interpolate='spline', bias=1, space="rgb")

	plot(grid.c, col=col.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) mtext("colonization probability", cex=cex.title)
	mtext(spText, side=2, cex=cex.title)
	plot(grid.e, col=ext.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("extinction probability", cex=cex.title)
	plot(grid.lam, col=pres.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("Probability of presence (CE Model)", cex=cex.title)
	plot(grid.sdm, col=sdm.colors(100), xaxt='n', yaxt='n', axis.args=leg.args)
	plotbg()
	if(i == 1) 	mtext("Probability of presence (SDM)", cex=cex.title)

}
dev.off()

