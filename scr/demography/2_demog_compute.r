#!/usr/bin/env Rscript

# qsub script for this one:
# module load R
# qsub -q qfat256 -l walltime=72:00:00 -l nodes=1:ppn=1 scr/demography/2_demog_compute.r

# computes number of species recruiting or dying by year in each plot
# requires:
#  	dat/demog/trees_wide.rds
#  	dat/demog/sample_years.rds

# produces:
#  	res/demography/demography.rds

library(foreach)
library(doParallel)
library(iterators)

numNodes = 48

outDir = file.path('res', 'demography')
suppressWarnings(dir.create(outDir, recursive=TRUE))

trees.wide = readRDS("dat/demog/trees_wide.rds")
sampleYears = readRDS("dat/demog/sample_years.rds")
speciesList = unique(trees.wide$id_spe)

arg = commandArgs(TRUE)
if(length(arg) > 0)
{
	sps = which(speciesList %in% arg)
	if(length(sps) == 0) 
	{
		warning("No species specified on command line; falling back to default")
	} else 
		speciesList = speciesList[which(speciesList %in% arg)]
}


trees.wide = trees.wide[trees.wide[,'id_spe'] %in% speciesList,]

# find all combinations of species and plot
spPlots = unique(trees.wide[,c('id_spe', 'plot_id')])


cl <- makeCluster(numNodes)
registerDoParallel(cl)
demog = foreach(spp = iter(spPlots, by='row'), .combine=rbind, .packages='foreach',
		.errorhandling = 'stop') %dopar%
{
	sp = spp[['id_spe']]
	pl = as.character(spp[['plot_id']])
	dat = with(trees.wide, trees.wide[id_spe == sp & plot_id == pl,])
	yrs = colnames(sampleYears)[which(sampleYears[pl,] > 0)]
	
	# little utility functions just to make the code below clearer:
	alive = function(plot.num, year)
		dat[dat$plot_id == plot.num,year] > 0
	not.alive = function(plot.num, year)
		dat[dat$plot_id == plot.num,year] == 0
		
	foreach(i=2:length(yrs), .combine=rbind, .errorhandling = 'remove') %do% {
		yr = yrs[i]
		yr_prev = yrs[i-1]

		# interval: the time interval since the previous time step
		# type: the rde type
		# recruit: the number of trees recruiting into the current time step since the previous
		# N.recruit: the sample size for recruits (i.e., the total alive in the current time step)
		# died: the number of trees that were alive in the previous time step and are dead in this one
		# N.died: the sample size for died; the number of trees that were alive in the previous time step
		data.frame(plot_id = as.integer(pl), species=sp, 
			interval=as.numeric(yr) - as.numeric(yr_prev), year = as.numeric(yr), 
			type=unique(dat[dat$plot_id == pl,'type']),		# should be length 1
			recruit = sum(alive(pl, yr) & not.alive(pl, yr_prev)),
			N.recruit = sum(alive(pl, yr)),
			died = sum(not.alive(pl, yr) & alive(pl, yr_prev)),
			N.died = sum(alive(pl,yr_prev)))
	}
}

## get rid of any rows where there was no recruitement/mortality and no trees
demog = with(demog, demog[recruit > 0 | died > 0,])

num = floor(runif(1,0,9999999))
#save the result
saveRDS(demog, file.path(outDir, paste0('demography',num, '.rds')))
stopCluster(cl)

