#!/usr/bin/env Rscript
library(lme4)
## library(rstan)
## rstan_options(auto_write = TRUE)
## options(mc.cores = parallel::detectCores())

## setwd("~/Dropbox/work/projects/STModel-Two-State_git")
species = read.table("dat/raw/dbh_trees_20151026.csv", header=TRUE, sep=';', dec='.', stringsAsFactors=FALSE)
species = species[complete.cases(species),]

# dbh filter to follow the original data
species = species[species$dbh > 127 & species$dbh < 9999,]


speciesList = unique(species$id_spe)
## speciesInfo = read.csv('dat/speciesInfo.csv')
source('scr/stm_functions.r')
rdeCols = c('#1f78b4', '#fb9a99')

# make a list of all unique points
plots = unique(species[,c('plot_id', 'longitude', 'latitude')])
coordinates(plots) = c('longitude', 'latitude')
proj4string(plots) = P4S.latlon

spIn = list()
for(spName in speciesList)
{
## 	info = speciesInfo[speciesInfo$spName == spName,]
## 	plLab = bquote(italic(.(as.character(info$genus))~.(as.character(info$species))))

	spGrid = readRDS(file.path('res','maps',paste0(spName,'_maps.rds')))
	
	# rde is fucked up for some reason; fix it here temporarily
	spGrid = within(spGrid,
	{
		rde[rde.present >= rde.contract & rde.present >= rde.expand ] = 0
		rde[rde.expand >= rde.contract & rde.expand > rde.present ] = 1
		rde[rde.contract > rde.present & rde.contract > rde.expand ] = 2
		rde[sdm < 0.1 & stm < 0.1 ] = NA
		
	})
	
	# make a raster of spGrid$rde
	rdeRas = make_raster(spGrid$rde, spGrid[,1:2])

	# extract the raster values at the points
	spPlot = data.frame(plot_id = plots$plot_id,
		type = extract(rdeRas, plots),
		id_spe = spName)
	spPlot = spPlot[complete.cases(spPlot),]
## 	spPlot = spPlot[spPlot$type != 1,]
	spIn[[spName]] = spPlot
}
spIn.df = do.call(rbind, spIn)

# set the in_range column in species to equal the point value from the extraction
#    might need a merge based on plot_id?
species.merged = merge(species, spIn.df, by=c('plot_id', 'id_spe'), all.x=TRUE)
species2 = species.merged
# drop type 1 (expansion)
species2 = species2[species2$type != 1,]
# keep only live trees
species2 = species2[species2$is_dead == 'f',]
## species2$is_dead[species2$is_dead == ""] = NA
## species2 = species2[complete.cases(species2),]
# print warning if not all cases are complete
if(sum(complete.cases(species2)) != nrow(species2))
{
	warning(paste(nrow(species2) - sum(complete.cases(species2)), "NA's found in species2 out of", nrow(species2), "rows"))
	species2 = species2[complete.cases(species2),]
}
species2$plot_id = factor(species2$plot_id)
species2$id_spe = factor(species2$id_spe)
species2$type = factor(species2$type)
species2$year_measured = factor(species2$year_measured)
## species2$is_dead = factor(species2$is_dead)


## par(las=2, cex.axis=0.5)
## box.locs = 1:20 + rep(0:9, each=2)
## boxplot(dbh ~ type + id_spe, data=species2,  
## 		notch=TRUE, range=1.5,outline=FALSE, col=rdeCols, boxwex=0.5, at=box.locs)
## 
## species2.subsample=species2.live[sample(1:nrow(species2.live),as.integer(nrow(species2.live)/10)),]	
## mod1 = lmer(dbh ~ type + (1|id_spe), data=species2.subsample)
## mod2 = lmer(dbh ~  type + (type|id_spe), data=species2.subsample)
## 
# fit mod2 in stan

## species2.stan = species2
## stdat = list(
## 	num_data_points = nrow(species2.stan),
## 	num_species = length(levels(species2.stan$id_spe)),
## 	dbh = species2.stan$dbh,
## 	species = as.integer(species2.stan$id_spe),
## 	type = as.integer(species2.stan$type) - 1)
## stdat2 = list(
## 	num_data_points = nrow(species2.stan),
## 	num_species = length(levels(species2.stan$id_spe)),
## 	num_years = length(levels(species2.stan$year_measured)),
## 	num_plots = length(levels(species2.stan$plot_id)),
## 	dbh = species2.stan$dbh,
## 	species = as.integer(species2.stan$id_spe),
## 	year = as.integer(species2.stan$year_measured),
## 	plot = as.integer(species2.stan$plot_id),
## 	type = as.integer(species2.stan$type) - 1)
## 
## cat('starting first stan model\n')
## stanMod = stan(file='scr/size_dist.stan', dat=stdat, iter=5000, chains=4)
## saveRDS(stanMod, 'res/dbhStanMod.rds')
## 
## cat('starting second stan model\n')
## stanMod2 = stan(file='scr/size_dist2.stan', dat=stdat2, iter=5000, chains=4)
## saveRDS(stanMod2, 'res/dbhStanMod2.rds')

gap=0.3
boxwidth = 0.4
between=5
nspecies=10
xlocs = rep(seq(0, by=between, length.out=nspecies), each=2) + rep(c(0,gap), nspecies)
xlim = range(xlocs) + 0.5*c(-between, between)
cex.pt = 0.7
rdeCols = c('#1f78b4', '#fb9a99')
cols = rep(rdeCols, nspecies)
cols=paste0(cols, "88")

## stan models were no good, too much data - try to get CIs via a bootstrap
# model 1: no random effects
m1 = lm(dbh ~ type*id_spe, data=species2)
nd = expand.grid(levels(species2$type), levels(species2$id_spe))
names(nd) = c('type', 'id_spe')
m1.pr = cbind(nd, predict(m1, newdata=nd, interval='prediction'))

par(mfrow=c(1,2), mar=c(0.5, 4, 0.5, 0.5))
plot(0,0, xlim=xlim, ylim=c(min(m1.pr$lwr), max(m1.pr$upr)), type='n', xaxt='n',xlab='',  ylab="dbh (mm)")
for(i in 1:nrow(m1.pr))
{
	xx = xlocs[i] + c(-boxwidth,boxwidth, boxwidth, -boxwidth)
	yy = rep(c(m1.pr$lwr[i], m1.pr$upr[i]), each=2)
	polygon(xx, yy, border=NA, col=cols[i])
	points(xlocs[i], m1.pr$fit[i], pch=16, cex=cex.pt)
}

inds = c(2, 12:20)
effs = c(coefficients(m1)[inds[1]], coefficients(m1)[inds[-1]] + coefficients(m1)[inds[1]])
cis = confint(m1, inds)
cis[2:nrow(cis),] = cis[2:nrow(cis),] + effs[1]

xl2 = seq(0, by=between, length.out = nspecies)
xlim2 = range(xl2) + 0.5*c(-between, between)
plot(0,0, xlim=xlim, ylim=c(min(cis[,1]), max(cis[,2])), type='n', xaxt='n',xlab='',  ylab=bquote(Delta*dbh~(mm)))
for(i in 1:nrow(cis))
{
	xx = xl2[i] + c(-boxwidth,boxwidth, boxwidth, -boxwidth)
	yy = rep(c(cis[i,1], cis[i,2]), each=2)
	polygon(xx, yy, border=NA, col=rdeCols[2])
	points(xl2[i], effs[i], pch=16, cex=cex.pt)
}



# model 2: random plot effect
m2 = lmer(dbh ~ type*id_spe + (1|plot_id), data=species2)

# model 3: year nested in plot
m3 = lmer(dbh ~ type*id_spe + (1|plot_id/year_measured), data=species2)

