#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

arg = commandArgs(TRUE)
modType = if('-m' %in% arg) 'mortality' else 'recruitment'

outDir = file.path('res', 'demography')
demog = readRDS(file.path(outDir, 'demography.rds'))
demog = demog[demog$interval <= 15,]
## demog = within(demog, {
## 	r.rate = recruit/N.recruit
## 	d.rate = died/N.died
## })

## bp.at = 1:20 + c(0, -0.5)
## par(mar=c(7,2,0.5,0.5))
## boxplot(r.rate ~ type+species, data=demog, las=3, cex.axis=0.5, boxwex = 0.4, at=bp.at)

## boxplot(d.rate ~ type+species, data=demog, las=3, cex.axis=0.5, boxwex = 0.4, at=bp.at)

## m1 = glm(cbind(recruit, N.recruit) ~ species*type, data=demog, family=binomial)
## library(lme4)
## m2 = glmer(cbind(recruit, N.recruit) ~ type + (type | species), data=demog, family=binomial) 
## m3 = glmer(cbind(recruit, N.recruit) ~ type + (1 | species), data=demog, family=binomial) 
## m4 = glmer(cbind(recruit, N.recruit) ~ type + (type | plot_id/species), data=demog, family=binomial) 
## 
## 
## m5 = glmer(cbind(recruit, N.recruit) ~ type + (1 | plot_id), data=demog[demog$species=='18032-ABI-BAL',], family=binomial) 
## m6 = glmer(cbind(recruit, N.recruit) ~ type + (type | plot_id), data=demog[demog$species=='18032-ABI-BAL',], family=binomial) 

stnResultPath.r = file.path('res', 'recruit_stan.rds')
stdat.r = with(demog[demog$N.recruit > 0,], list(
	num_data_points = length(recruit),
	num_species = length(levels(species)),
	recruit = recruit,
	N = N.recruit,
	species = as.integer(species),
	type = as.integer(type) - 1,
	interval = interval))

stnResultPath.m = file.path('res', 'recruit_stan.rds')
stdat.m = with(demog[demog$N.died > 0,], list(
	num_data_points = length(died),
	num_species = length(levels(species)),
	recruit = died,
	N = N.died,
	species = as.integer(species),
	type = as.integer(type) - 1,
	interval = interval))

results = list(data.orig = demog,
		species = data.frame(name = unique(demog$species), 
				value = unique(as.integer(demog$species))),
		type = data.frame(name = unique(demog$type), 
				value = unique(as.integer(demog$type)-1)))

if(modType == 'recruitment')
{
	cat('launching stan recruitment model\n')
	stMod = stan(file='scr/demography/demography.stan', dat=stdat.r, iter=5000, chains=4)
	results$data.stan = stdat.r
	fname = 'res/demography/recruitment_stan.rds'
} else {
	cat('launching stan mortality model\n')
	stMod = stan(file='scr/demography/demography.stan', dat=stdat.m, iter=5000, chains=4)
	results$data.stan = stdat.m
	fname = 'res/demography/mortality_stan.rds'
}
results$results = stMod
saveRDS(results, fname)
