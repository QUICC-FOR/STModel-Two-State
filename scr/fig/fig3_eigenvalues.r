#!/usr/bin/env Rscript

## figures produced:
##    img/figs/fig3.png

library(numDeriv)
## library(ellipse)
speciesInfo = read.csv('dat/speciesInfo.csv')
climScale = readRDS('dat/climate_scaling.rds')
dldt.path = file.path('res','dldt.rds')

north.col = '#4575b4'
south.col = '#d73027'
lam_t = function(T, pars, P)
{
	Cl = pars[1] + pars[2]*T + pars[3]*P + pars[4]*T^2 + pars[5]*P^2
	El = pars[6] + pars[7]*T + pars[8]*P + pars[9]*T^2 + pars[10]*P^2
	plogis(Cl) - plogis(El)
}


get_derivs = function(curPars, temp, precip, stepSize = 0.25)
{
	starts = seq(min(temp), max(temp), stepSize)

	roots = sapply(starts, function(st) 
		{
			tryCatch(
			uniroot(lam_t, pars = curPars, P = precip, lower = st, upper = st + stepSize)$root,
			error=function(e) NA,
			warning=function(w) NA)
		})
	roots = roots[!is.na(roots)]
	roots = unique(roots)

	dldt = grad(lam_t, roots, pars = curPars, P = precip)
	type = c(1,2)
	if(length(roots) == 1) type = -1

	matrix(c(dldt, roots, type), nrow=length(roots), ncol=3)
}

pct.done = function(pct, overwrite = TRUE, pad='', digits = 0)
{
	if(overwrite) cat('\r')
	cat(paste0(pad, round(pct, digits), '%'))
	flush.console()
}

compute_species_dldt = function(spName)
{
	# get the calibration range - the range limits must be within this range
	calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds')))
	precip = median(calibDat$tot_annual_pp)

	# get posterior
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_posterior.rds'))
	spPosterior = readRDS(posteriorFname)[['0']]

	outputSteps = seq(floor(0.01*nrow(spPosterior)), nrow(spPosterior), length.out=100)

	result = matrix(NA, nrow = (2*nrow(spPosterior)+4), ncol=3)
	curRow = 1
	cat(spName,'\n')
	pct.done(0, FALSE, "computing dldt: ")
	for(i in 1:nrow(spPosterior))
	{
		curPars = spPosterior[i,]
		dldt = get_derivs(curPars, calibDat$annual_mean_temp, precip)
		result[curRow:(curRow + nrow(dldt) - 1),] = dldt
		curRow = curRow + nrow(dldt)
		if(curRow > nrow(result)) stop("Something is wrong")
		if(i %in% outputSteps)
			pct.done(100* i / nrow(spPosterior), pad = 'computing dldt: ')
	}
	cat('\n')
	result = result[1:(curRow - 1),]

	# classify single roots as belonging either upper or lower range limit
	singles = which(result[,3]==-1)
	lowers = which(result[,3]==1)
	uppers = which(result[,3]==2)
	lowdist = abs(result[singles,2] - mean(result[lowers,2]))
	updist = abs(result[singles,2] - mean(result[uppers,2]))
	result[singles,3][lowdist <= updist] = 1
	result[singles,3][lowdist > updist] = 2

	# recompute after reclassifying
	lowers = which(result[,3]==1)
	uppers = which(result[,3]==2)
	result[,1] = abs(result[,1])

	lapply(list(north=lowers,south=uppers), function(ind) 
	{
		c(T=mean(result[ind,2]), dldt=mean(result[ind,1]),
		 cor=cor(result[ind,2],result[ind,1]), sd.T = sd(result[ind,2]), 
		 sd.dldt=sd(result[ind,1]))
	})
}

if(!file.exists(dldt.path))
{
	spList = as.list(speciesList)
	names(spList) = speciesList
	dldt = lapply(spList, compute_species_dldt)

	dldt.df = data.frame(spcode=character(0), genus=character(0), species=character(0),
			boundary=character(0), T=numeric(0), dldt=numeric(0), cor=numeric(0), 
			sd.T = numeric(0), sd.dldt=numeric(0))
	for(sp in speciesList)
	{
		info = speciesInfo[speciesInfo$spName == sp,]
		for(bound in c("north", "south"))
		{
			X = dldt[[sp]][[bound]]
			dldt.df = rbind(dldt.df,data.frame(spcode=sp, genus=info$genus, species=info$species,
					boundary=bound, T=X[1], dldt=X[2], cor=X[3], sd.T=X[4], sd.dldt=X[5]))
		}
	}
	rownames(dldt.df) = 1:nrow(dldt.df)
	# fix picmar - the boundary gets incorrectly classified
	dldt.df = dldt.df[-which(dldt.df$spcode == '183302-PIC-MAR' & dldt.df$boundary=='south'),]
	dldt.df$boundary[dldt.df$spcode == '183302-PIC-MAR'] = 'south'

	# transform the x-coordinates back to the original range
	dldt.df = within(dldt.df, {
		T = (T * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp']})
	## 	T.lower = (T.lower * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp']
	## 	T.upper = (T.upper * climScale$scale['annual_mean_temp']) + climScale$center['annual_mean_temp']})

	dldt.df$color = north.col
	dldt.df$color[dldt.df$boundary == 'south'] = south.col

	saveRDS(dldt.df,dldt.path)
} else
{
	dldt.df = readRDS(dldt.path)
}


# make the actual figure

dpi = 600
figure.width = 4
figure.height = 4.5
filename = file.path('img', 'figs', 'fig3.png')
fontsize=12
png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
	file=filename, pointsize=fontsize, res=dpi)
par(mar=c(3,3,0.5,1.5), mgp=c(1.5,0.5,0), tcl=-0.2, cex.axis=0.6, cex.lab=0.7, xpd=NA)

plot(0,0,type='n', xlab="Mean Annual Temperature (°C)", ylab=expression("|"*partialdiff*lambda/partialdiff*T*"|"),
		xlim=c(-5,25), ylim=c(0,0.6), bty='n')
for(sp in dldt.df$spcode)
{
	info = speciesInfo[speciesInfo$spName == sp,]
	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))
	xs = dldt.df[dldt.df$spcode == sp & dldt.df$boundary == 'south',]
	xn = dldt.df[dldt.df$spcode == sp & dldt.df$boundary == 'north',]
	
	# the ellipses show multivariate 95% confidence limits
## 	el = ellipse(xs$cor, scale=c(xs$sd.T, xs$sd.dldt), centre=c(xs$T, xs$dldt))
## 	polygon(el[,1], el[,2], bg="NA", col=paste0(xs$color,'22'), border=NA)
	if(sp != '183302-PIC-MAR')
	{
		segments(xn$T, xn$dldt, xs$T, xs$dldt, lty=1, lwd=0.7, col="#555555")
## 		el = ellipse(xn$cor, scale=c(xn$sd.T, xn$sd.dldt), centre=c(xn$T, xn$dldt))
## 		polygon(el[,1], el[,2], bg="NA", col=paste0(xn$color,'22'), border=NA)
	}
	a = text(xs$T, xs$dldt, plLab, pos=4, cex=0.55, col="#444444")
}
with(dldt.df, points(T, dldt, col=color, pch=16))
## legend(-4.5, 0.6, legend=c("Northern range limit", 
## 		"Southern range limit"), col=c(north.col, south.col), 
## 		text.col = c(north.col, south.col), pch=16, bty='n', cex=0.7)
text(-6, -0.01, expression(bold(Northern~range~limits)), col=north.col, cex=0.6, pos=4)
text(8, -0.01, expression(bold(Southern~range~limits)), col=south.col, cex=0.6, pos=4)
arrows(-10,0.45, -10, 0.6, length=.07, lwd=1.5)
text(-11.5,0.45, "Faster dynamics", srt=90, cex=0.5, pos=4)
arrows(-10,0.15, -10, 0, length=.07, lwd=1.5)
text(-9.5,0.15, "Slower dynamics", srt=90, cex=0.5, pos=2)
