#!/usr/bin/Rscript

library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
argList = parser$parse_args()
spName = argList$species
# set seed - drawn from sample(1:1e6, 1)
set.seed(588533)

#-------------------
#  Load and prepare data
#-------------------

infile = paste("dat/transition_twostate_", spName, ".rdata", sep="")
load(infile)

# get and scale the climate variables
colnames(stateData)[which(colnames(stateData) == spName)] = "presence"
climVars = which(!(colnames(stateData) %in% c('presence', 'plot_id', 'year_measured', 'lat', 'lon')))
climVarNames = colnames(stateData)[climVars]
stateData.clim.scaled = scale(stateData[,climVars])

# save the transformation as a function
scalefunc = function(sc_mean, sc_sd) {
	return(list(
		scale = function(dat)
		{
			for(nm in names(sc_mean)) dat[,nm] = (dat[,nm] - sc_mean[nm])/sc_sd[nm]
			return(dat)
		},
		unscale = function(dat)
		{
			for(nm in names(sc_mean)) dat[,nm] = dat[,nm]*sc_sd[nm] + sc_mean[nm]
			return(dat)
		},
		means = sc_mean,
		sds = sc_sd))
}
clim.scale = scalefunc(attr(stateData.clim.scaled, "scaled:center"), attr(stateData.clim.scaled, "scaled:scale"))

stateData.scaled = cbind(stateData[,-climVars], stateData.clim.scaled)
transitionData.scaled = clim.scale$scale(transitionData)

# subset the data for the SDM
# select only a single row for each plot (to avoid too much spatial duplication)
rows = sapply(unique(stateData.scaled$plot_id), function(i) {
	candidates = which(stateData.scaled$plot_id == i)
	if(length(candidates) == 1) candidates else sample(candidates, 1)
	})
stateData.subset = stateData.scaled[rows,]


save(stateData.subset, transitionData.scaled, spName, climVars, 
		climVarNames, clim.scale, file=paste("dat/", spName, "_processed.rdata", sep=""))


climGrid.raw = read.csv("dat/SDMClimate_grid.csv")
climGrid.raw = climGrid.raw[complete.cases(climGrid.raw),]
# apply transformation to the climate grid
climGrid.scaled = clim.scale$scale(climGrid.raw)
saveRDS(climGrid.scaled, file = paste("dat/", spName, "_climData_scaled.rds", sep=""))

climGrid.past = read.csv("dat/2states_past_climate_1901-1930.csv")
climGrid.past = climGrid.past[complete.cases(climGrid.past),]
climGrid.past.scaled = clim.scale$scale(climGrid.past)
saveRDS(climGrid.past.scaled, file = "dat/climPast_scaled.rds")
