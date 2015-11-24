#!/usr/bin/env Rscript

## figures produced:
##    img/figs/fig3.png

library(numDeriv)
library(rootSolve)
library(foreach)
library(iterators)
library(coda)
library(doParallel)
registerDoParallel(cores=detectCores())

speciesInfo = read.csv('dat/speciesInfo.csv')
climScale = readRDS('dat/clim/climate_scaling.rds')
dldt.path = file.path('res','dldt.rds')
speciesList = readRDS('dat/speciesList.rds')
posterior.n = 1000
climGrid = readRDS(file.path('dat', 'clim', 'climateGrid_scaled.rds'))

recompute = FALSE
arg = commandArgs(trailingOnly = TRUE)
if('--recompute' %in% arg | '-r' %in% arg) recompute = TRUE

north.col = '#4575b4'
south.col = '#d73027'
lam_t = function(T, pars, P)
{
	Cl = pars[1] + pars[2]*T + pars[3]*P + pars[4]*T^2 + pars[5]*P^2
	El = pars[6] + pars[7]*T + pars[8]*P + pars[9]*T^2 + pars[10]*P^2
	plogis(Cl) - plogis(El)
}


get_derivs = function(curPars, temp.range, precip, stepSize = 0.25)
{
	roots = uniroot.all(lam_t, pars=curPars, P = precip, lower=temp.range[1], upper=temp.range[2])
	if(length(roots) > 0)
	{
		dldt = grad(lam_t, roots, pars = curPars, P = precip)
		type = sapply(roots, function(r)
			# if we get colder (root - 0.01) and lambda is negative, then it is a northern limit
			ifelse(lam_t(r-0.01, curPars, precip) < 0, "northern", "southern"))
	} else {
		dldt = roots = type = NA
	}

	data.frame(dldt=dldt, root=roots, boundary=type, stringsAsFactors = FALSE)
}

## pct.done = function(pct, overwrite = TRUE, pad='', digits = 0)
## {
## 	if(overwrite) cat('\r')
## 	cat(paste0(pad, round(pct, digits), '%'))
## 	flush.console()
## }

compute_species_dldt = function(spName)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	# get the median precipitation experienced at the range boundary
	# we use the 'uncertainty range' where p is between 0.05 and 0.95
	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	spGrid = merge(spGrid, climGrid[,c('x','y', 'annual_mean_temp', 'tot_annual_pp')])
	precip = median(spGrid$tot_annual_pp[spGrid$stm>0.05 & spGrid$stm < 0.95], na.rm=TRUE)
	trange = range(spGrid$annual_mean_temp)
	
	# get posterior
	posteriorFname = file.path('res', 'posterior', paste0(spName, '_0_samples.rds'))
	samples = readRDS(posteriorFname)
	samples = samples[seq(1, nrow(samples), length.out=posterior.n),]

	dldt = foreach(curPars=iter(samples, by='row'), .combine=rbind) %do%
			get_derivs(curPars, trange, precip)
	dldt$dldt = abs(dldt$dldt)
	dldt$precip = precip
	dldt$spName = spName
	dldt$type = info$type


	# picmar and pinban boundaries are misclassified for some reason
	# drop the fake northern boundary
	if(spName == "183319-PIN-BAN" | spName == "183302-PIC-MAR")
		# drop all fake northern boundaries and relabel the northern boundary
		dldt = dldt[-which(dldt$boundary == 'northern'),]


	dldt[complete.cases(dldt),]
## 	if(nrow(dldt) > 0) {
## 		data.frame(spcode = spName, genus=info$genus, species=info$species,
## 			boundary = unique(dldt$type),
## 			dldt.mean=tapply(dldt$dldt, dldt$type, mean),
## 			dldt.sd = tapply(dldt$dldt, dldt$type, sd),
## 			temp.mean=tapply(dldt$root, dldt$type, mean),
## 			temp.sd=tapply(dldt$root, dldt$type, sd),
## 			dldt.cor = cor(dldt$root, dldt$dldt))
## 	} else {
## 		data.frame(spcode = spName, genus=info$genus, species=info$species, boundary = NA, 
## 			dldt.mean=NA, dldt.sd = NA, temp.mean=NA, temp.sd=NA, dldt.cor = NA)
## 	}
	
}

if(!file.exists(dldt.path) | recompute)
{
	dldt = foreach(spName = speciesList, .combine = rbind) %dopar%
		 compute_species_dldt(spName)
## 	rownames(dldt) = 1:nrow(dldt)

	# transform the x-coordinates back to the original range
	dldt = within(dldt, { root = (root * climScale$scale['annual_mean_temp']) + 
				climScale$center['annual_mean_temp']})

	saveRDS(dldt,dldt.path)
} else
{
	dldt = readRDS(dldt.path)
}




# if we go down the lmer path (or stan) this is the model I liked the best
## drp = which(dldt$spName == "183319-PIN-BAN" | dldt$spName == "183302-PIC-MAR")
## will have to simulate northern range boundaries for larlar pinban and picmar
## dldt$boundary = factor(dldt$boundary)
## dldt$spName = factor(dldt$spName)
## summary(lmer(dldt ~ boundary + (1+boundary|spName), data=dldt[-drp,]))

# make the actual figure

dpi = 600
figure.width = 5.5/2.54
figure.height = 6/2.54
filename = file.path('img', 'figs', 'fig3.pdf')
fontsize=9
pdf(width=figure.width, height=figure.height, file=filename, pointsize=fontsize)
## png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
## 	file=filename, pointsize=fontsize, res=dpi)
par(mar=c(3,3,0.5,1.5), mgp=c(1.5,0.5,0), tcl=-0.2, cex.axis=0.7, cex.lab=0.8, xpd=NA)



dldt.fig = aggregate(.~boundary+spName+type, data=dldt, FUN=mean)
# exclude species that only have one range boundary
exclusions = c("183319-PIN-BAN", "183412-LAR-LAR", "183302-PIC-MAR")
dldt.fig = dldt.fig[!(dldt.fig$spName %in% exclusions),]

lab.offsets = data.frame(spCode=unique(dldt.fig$spName), x=0, y=0)
lab.offsets[lab.offsets$spCode == "28731-ACE-SAC", c('x', 'y')] = c(0.5, 0.0)
lab.offsets[lab.offsets$spCode == "183385-PIN-STR", c('x', 'y')] = c(-0.5, -0.02)
lab.offsets[lab.offsets$spCode == "19408-QUE-RUB", c('x', 'y')] = c(0.2, 0.04)
lab.offsets[lab.offsets$spCode == "19489-BET-PAP", c('x', 'y')] = c(-0.5, 0.02)
lab.offsets[lab.offsets$spCode == "18032-ABI-BAL", c('x', 'y')] = c(-0.5, 0.02)
lab.offsets[lab.offsets$spCode == "195773-POP-TRE", c('x', 'y')] = c(-.3, -0.02)
lab.offsets[lab.offsets$spCode == "183375-PIN-RES", c('x', 'y')] = c(-0.25, 0.01)
lab.offsets[lab.offsets$spCode == "19481-BET-ALL", c('x', 'y')] = c(-3.5, -0.035)
lab.offsets[lab.offsets$spCode == "183397-TSU-CAN", c('x', 'y')] = c(-0.5, 0.015)
lab.offsets[lab.offsets$spCode == "183295-PIC-GLA", c('x', 'y')] = c(-0.5, -0.005)
lab.offsets[lab.offsets$spCode == "28728-ACE-RUB", c('x', 'y')] = c(-0.5, -0.01)
lab.offsets[lab.offsets$spCode == "19462-FAG-GRA", c('x', 'y')] = c(-0.5, -0.01)


pch.types = c(temperate=16, boreal=17, transitional=18)
cex.pt.types = c(temperate=0.8, boreal=0.8, transitional=0.8)

yl = c(0,1)
plot(0,0,type='n', xlab="Mean Annual Temperature (°C)", ylab=expression(partialdiff*lambda/partialdiff*T),
		xlim=c(-5,20), ylim = yl, bty='n')
for(sp in unique(as.character(dldt.fig$spName)))
{
	dat = dldt.fig[dldt.fig$spName == sp,]
	info = speciesInfo[speciesInfo$spName == sp,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus), 1, 1), '.'))~.(as.character(info$species))))
	pch = pch.types[as.character(info$type)]
	cex.pt = cex.pt.types[as.character(info$type)]
	oset.x = lab.offsets[lab.offsets$spCode == sp, 'x']
	oset.y = lab.offsets[lab.offsets$spCode == sp, 'y']
	
	xs = dat[dat$boundary == 'southern',]
	xn = dat[dat$boundary == 'northern',]
	tl.x = xs$root
	tl.y = xs$dldt
	segments(xn$root, xn$dldt, xs$root, xs$dldt, lty=1, lwd=0.7, col="#555555")
	points(xs$root, xs$dldt, col=south.col, pch=pch, cex=cex.pt)
	points(xn$root, xn$dldt, col=north.col, pch=pch, cex=cex.pt)
	text(tl.x+oset.x, tl.y+oset.y, plLab, pos=4, cex=0.45, col="#444444")
}
legend(-4.5, 0.8, legend=c("Temperate", "Transitional", "Boreal"), 
		pch=c(pch.types[c('temperate', 'transitional', 'boreal')]), 
		col='black', bty='n', cex=0.6)
## text(-3, 0, expression(bold(Northern~range~limits)), col=north.col, cex=0.6, pos=4)
## text(7, 0, expression(bold(Southern~range~limits)), col=south.col, cex=0.6, pos=4)
par(xpd=NA)
arrow.x = -11
atext.oset=-2.5
mid = 0.5
a.len = 0.25
a.oset = 0.2
slow.text.yoset = 0.3
arrows(arrow.x,mid+a.oset, arrow.x, mid+a.oset+a.len, length=.05, lwd=1.25)
text(arrow.x+atext.oset,mid+a.oset, "Faster dynamics", srt=90, cex=0.5, pos=4)
arrows(arrow.x,mid-a.oset, arrow.x, mid-a.oset-a.len, length=.05, lwd=1.25)
text(arrow.x+atext.oset,mid-a.oset-slow.text.yoset, "Slower dynamics", srt=90, cex=0.5, pos=4)
dev.off()
