## figures produced:
##    img/resp_curves.png (all species)

library(coda)
speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv')

source('scr/stm_functions.r')

# colors
col.color = '#377eb8'
col.bgcolor = '#377eb855'
ext.color = '#e41a1c'
ext.bgcolor = '#e41a1c55'


dpi = 600
figure.width = 6.5
figure.height = 8
filename = file.path('img', 'resp_curves.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mfrow=c(5, 4), bty='n', mar=c(2,2,0,0), mgp=c(1,0.25,0), oma=c(2,2,0,0), tcl=-0.2, cex.axis=0.5)
cex.xtitles = 0.5
cex.ytitle = 0.75
cex.title = 0.5
line.title = -1
rnum = 1
cnum = 1
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	rc = readRDS(file.path('res','resp_curve', paste0(spName, '_respCurve.rds')))
	plot(rc$temp, rc$col.temp, col=col.color, type='l', xlab='', ylab='', ylim=c(0,info$rc_tymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$temp, rev(rc$temp)), c(rc$col.temp.lower, rev(rc$col.temp.upper)), col=col.bgcolor, border=NA)
	lines(rc$temp, rc$ext.temp, col=ext.color) 
	polygon(c(rc$temp, rev(rc$temp)), c(rc$ext.temp.lower, rev(rc$ext.temp.upper)), col=ext.bgcolor, border=NA)
	if(rnum == 5)
		mtext("Mean Annual Temperature (°)", side = 1, cex=cex.xtitles, line=2)
	plot(rc$precip, rc$col.precip, col=col.color, type='l', xlab='', ylab='', ylim=c(0, info$rc_pymax))
	mtext(plLab, side=3, cex=cex.title, line=line.title)
	polygon(c(rc$precip, rev(rc$precip)), c(rc$col.precip.lower, rev(rc$col.precip.upper)), col=col.bgcolor, border=NA)
	lines(rc$precip, rc$ext.precip, col=ext.color) 
	polygon(c(rc$precip, rev(rc$precip)), c(rc$ext.precip.lower, rev(rc$ext.precip.upper)), col=ext.bgcolor, border=NA)
	if(rnum == 5)
		mtext("Mean Annual Precipitation (mm)", side = 1, cex=cex.xtitles, line=2)
	if(cnum == 2)
	{
		cnum = 1
		rnum = rnum+1
	} else
		cnum = cnum + 1
}
mtext(expression("Probability of " * phantom(bold(" colonization")) * "/" * 
		phantom(bold("extinction"))),col="black", side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * bold(" colonization") * phantom("/") * 
		phantom(bold("extinction"))),col=col.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
mtext(expression(phantom("Probability of ") * phantom(bold(" colonization")) * phantom("/") * 
		bold("extinction")),col=ext.color, side=2, outer=TRUE, cex=cex.ytitle, line=0)
dev.off()



stmCols = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
sdmCols = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
	interpolate='spline', bias=2, space="rgb")
stmVarCols = colorRampPalette(c("#ffffff", "#bdc9e1", "#045a8d", "#33338d", "#ffff88"), 
		interpolate='spline', bias=2, space="rgb")
rdeCols = c('#1f78b4', '#b2df8a', '#fb9a99')

spGrid = readRDS(file.path('res','maps',paste0(spName,'_maps.rds')))
## plot_sdm = function(sdmDat, coords, sdm.col, legend=FALSE, add=FALSE, plot.ocean=TRUE, ...)
quartz(w=10,h=3)
par(mfrow=c(1,4), mar=c(1,1,0,0), oma=c(0,0,1,1))
plot_sdm(spGrid$stm, spGrid[,1:2], stmCols(100))
plot_sdm(spGrid$stm.var, spGrid[,1:2], stmVarCols(100))
plot_sdm(spGrid$sdm, spGrid[,1:2], sdmCols(100))
plot_sdm(spGrid$rde, spGrid[,1:2], rdeCols)






## speciesList = c('28728-ACE-RUB', '28731-ACE-SAC','28728-ACE-RUB', '28731-ACE-SAC','28728-ACE-RUB', '28731-ACE-SAC','28728-ACE-RUB', '28731-ACE-SAC','28728-ACE-RUB', '28731-ACE-SAC')


for(spName in speciesList)
{


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
	
	
	cat.colors = c('#1f78b4', '#b2df8a', '#fb9a99')
	diseq = make_raster(spGrid$diseq, spGrid[,1:2], P4S.latlon, stmMapProjection)
	plot(diseq, col=cat.colors, xaxt='n', yaxt='n', legend=FALSE)
	plotbg()

dev.off()
		
		
	
}

