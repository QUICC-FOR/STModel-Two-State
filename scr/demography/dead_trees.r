#!/usr/bin/env Rscript
library(reshape2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('scr/stm_functions.r')

stnResultPath = file.path('res', 'deadStanMod.rds')
species = read.table("dat/raw/dbh_trees_20151026.csv", header=TRUE, sep=';', dec='.', stringsAsFactors=FALSE)
species = species[complete.cases(species),]

# dbh filter to follow the original data
species = species[species$dbh > 127 & species$dbh < 9999,]

speciesList = unique(species$id_spe)

# make a list of all unique points
plots = unique(species[,c('plot_id', 'longitude', 'latitude')])
coordinates(plots) = c('longitude', 'latitude')
proj4string(plots) = P4S.latlon

spIn = list()
for(spName in speciesList)
{
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
species.merged = merge(species, spIn.df, by=c('plot_id', 'id_spe'), all.x=TRUE)
species2 = species.merged
# drop type 1 (expansion)
species2 = species2[species2$type != 1,]
species2$is_dead[species2$is_dead == ''] = NA
if(sum(complete.cases(species2)) != nrow(species2))
{
	warning(paste(nrow(species2) - sum(complete.cases(species2)), "NA's found in species2 out of", nrow(species2), "rows"))
	species2 = species2[complete.cases(species2),]
}
species2$plot_id = factor(species2$plot_id)
species2$id_spe = factor(species2$id_spe)
species2$type = factor(species2$type)
species2$year_measured = factor(species2$year_measured)


# note to self; including random effects for plot and year improves the model, but it's not
# really necessary; the effects are the same no matter what
## m1 = glmer(cbind(t, t+f) ~ type + (1|id_spe), data=dead.noplots, family=binomial)
## m2 = glmer(cbind(t, t+f) ~ type + (type|id_spe), data=dead.noplots, family=binomial)

# fit m2 within stan

# now look at % dead by species in type 0 and type 1 plots
species2$weight = 1
temp = melt(species2, id.vars=c('id_spe', 'is_dead', 'type'), measure.vars='weight')
dead.noplots = dcast(temp, id_spe + type ~ is_dead, fun.aggregate=sum)

stdat = list(
	num_data_points = nrow(dead.noplots),
	num_species = length(levels(dead.noplots$id_spe)),
	dead = dead.noplots$t,
	alive = dead.noplots$f,
	species = as.integer(dead.noplots$id_spe),
	type = as.integer(dead.noplots$type) - 1)

stanMod = stan(file='scr/dead_trees.stan', dat=stdat, iter=5000, chains=3)
stnResults = list(
	data.orig = dead.noplots,
	data.stan = stdat,
	species = data.frame(name = unique(dead.noplots$id_spe), 
				value = unique(as.integer(dead.noplots$id_spe))),
	results = stanMod)

saveRDS(stnResults, 'res/deadStanMod.rds')


