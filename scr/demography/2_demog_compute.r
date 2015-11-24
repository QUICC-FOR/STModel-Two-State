#!/usr/bin/env Rscript

# computes number of species recruiting or dying by year in each plot
# requires:
#  	dat/demog/trees_wide.rds
#  	dat/demog/sample_years.rds

# produces:
#  	res/demography/demography.rds

library(foreach)
library(doParallel)

numNodes = 40

outDir = file.path('res', 'demography')
suppressWarnings(dir.create(outDir, recursive=TRUE))

trees.wide = readRDS("dat/demog/trees_wide.rds")
sampleYears = readRDS("dat/demog/sample_years.rds")
speciesList = unique(trees.wide$id_spe)


# find all combinations of species and plot
spPlots = unique(trees.wide[,c('id_spe', 'plot_id')])

cl <- makeCluster(numNodes)
registerDoParallel(cl)
demog = foreach(spp = iter(spPlots, by='row'), .combine=rbind, .packages='foreach') %dopar%
{
	sp = spp['id_spe']
	pl = spp['plot_id']
	dat = with(trees.wide, trees.wide[id_spe == sp & plot_id = pl,])
	yrs = colnames(sampleYears)[which(sampleYears[pl,] > 0)]
	res = foreach(i=2:length(yrs), .combine=rbind) %do% {
		yr = yrs[i]
		yr_prev = yrs[i-1]
		data.frame(plot_id = as.integer(pl), species=sp, 
			interval=as.numeric(yr) - as.numeric(yr_prev), year = as.numeric(yr), 
			type=unique(dat[dat$plot_id == pl,'type']),		# should be length 1
			recruit = sum(dat[dat$plot_id == pl,yr] > 0 & dat[dat$plot_id == pl,yr_prev] == 0),
			N.recruit = sum(dat[dat$plot_id == pl,yr] > 0),
			died = sum(dat[dat$plot_id == pl,yr] == 0 & dat[dat$plot_id == pl,yr_prev] > 0),
			N.died = sum(dat[dat$plot_id == pl,yr_prev] > 0))
	}
		# get rid of any rows where there was no recruitement/mortality and no trees
		res[!(res$recruit == 0 & res$died == 0),]
	}
}

#save the result
saveRDS(demog, file.path(outDir, 'demography.rds'))
stopCluster(cl)







