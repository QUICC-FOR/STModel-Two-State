#!/usr/bin/env Rscript

# does basic data reshaping to prepare for the demography analysis
# requires:
#   dat/raw/dbh_trees_20151026.rds
#         *or*
#   dat/raw/dbh_trees_20151026.csv

# creates:
#   dat/raw/dbh_trees_20151026.rds (if not present)
#  	dat/demog/trees_wide.rds
#  	dat/demog/sample_years.rds

library(rgdal)
library(foreach)
load("dat/map_projections.rdata")
source("scr/stm_functions.r")

library(reshape2)

suppressWarnings(dir.create(file.path('dat', 'demog'), recursive=TRUE))

cat("reading raw data\n")
if(!file.exists("dat/raw/dbh_trees_20151026.rds"))
{
	trees = read.table("dat/raw/dbh_trees_20151026.csv", header=TRUE, sep=';', dec='.', stringsAsFactors=FALSE)
	saveRDS(trees, "dat/raw/dbh_trees_20151026.rds")
} else {
	trees = readRDS("dat/raw/dbh_trees_20151026.rds")
}

speciesList = unique(trees$id_spe)

# properly assign NAs
trees$dbh[trees$dbh == 9999] = NA
trees$is_dead[trees$is_dead == ''] = NA
trees = trees[complete.cases(trees),]

# find plots that only have one sample
plFreqs = rowSums(with(trees, 1*(table(plot_id, year_measured)>0)))
singlePlots = as.numeric(names(plFreqs[which(plFreqs == 1)]))

# filters:
#   remove all NAs
#   plots with only one sample
#   dbh filter to follow the original data (>127 mm)
#   also drop dead trees
trees = trees[!(trees$plot_id %in% singlePlots) & trees$dbh > 127 & 
		trees$is_dead == 'f', !(colnames(trees) %in% c('is_dead'))]

# verify that trees are only ever observed in a single plot --- TRUE
## tr = trees$tree_id
## trPl = paste(trees$tree_id, trees$plot_id, sep='-')
## length(unique(tr)) == length(unique(trPl))





# figure out the category (present, expand, contract) for each trees/plot combo
# use the rde maps to do this, then merge back to the trees table
# warning: SLOW
# make a list of all unique points
modName = '0'
plots = unique(trees[,c('plot_id', 'longitude', 'latitude')])
coordinates(plots) = c('longitude', 'latitude')
proj4string(plots) = P4S.latlon

cat("computing plot categories\n")
spCategories = foreach(spName = speciesList, .combine=rbind, .final=data.frame) %do%
{
	spGrid = readRDS(file.path('res','maps',paste0(spName, '_', modName, '_maps.rds')))

	# make a raster of spGrid$rde
	rdeRas = make_raster(spGrid$rde, spGrid[,1:2])

	# extract the raster values at the points
	spPlot = data.frame(plot_id = plots$plot_id,
		type = extract(rdeRas, plots),
		id_spe = spName)
	spPlot = spPlot[complete.cases(spPlot),]
	spPlot
}

# we only care about present (category 0) and contraction (category 2)
spCategories = spCategories[spCategories$type != 1,]
spCategories$type = factor(spCategories$type, labels=c('present', 'contract'))

# merge the two - all is false because we only want rows where the category is present or contract
cat("merging categories\n")
trees = merge(trees, spCategories, by=c('plot_id', 'id_spe'), all.x=FALSE)


# cast into wide format
cat("casting data\n")
trees.wide = dcast(trees, tree_id + plot_id + id_spe + type ~ year_measured, value.var='dbh')

# set DBH to 0 for years where a tree was missing but the plot was sampled
# find years in which each plot was sampled and drop any with only one sample
cat("finding sampled years\n")
sampleYears = with(trees, table(trees$plot_id, trees$year_measured))
singlePlots = as.numeric(rownames(sampleYears)[which(rowSums(sampleYears > 0) < 2)])
trees.wide = trees.wide[!(trees.wide$plot_id %in% singlePlots),]

for(yr in colnames(sampleYears))
{
	pls = as.integer(rownames(sampleYears)[which(sampleYears[,yr] > 0)])
	tw.rows = which(trees.wide$plot_id %in% pls & is.na(trees.wide[,yr]))
	trees.wide[tw.rows, yr] = 0
}


saveRDS(trees.wide, "dat/demog/trees_wide.rds")
saveRDS(sampleYears, "dat/demog/sample_years.rds")

