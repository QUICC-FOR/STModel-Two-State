## library(foreach)
## library(doParallel)
library(reshape2)
## registerDoParallel(cores=2)
setwd("/Users/mtalluto/Dropbox/work/projects/STModel-Two-State_git/")

# process temp plot climate
climScale = readRDS("dat/climate_scaling.rds")
tp_clim = read.csv("dat/raw/temp_plots/tp_climData_reshaped.csv", dec=',', sep=';')
colsToKeep = c('plot_id', 'year_measured')
tp_clim = cbind(tp_clim[,colsToKeep], scale(tp_clim[,names(climScale$center)], 
		center=climScale$center, scale=climScale$scale))
tp_clim = tp_clim[complete.cases(tp_clim),]


tempPlots = read.csv("dat/raw/temp_plots/tp_plotinfoData.csv", dec='.', sep=';')

tp_species = read.csv("dat/raw/temp_plots/tp_treeData_allSpecies.csv", dec='.', sep=';')

treeDat = merge(tp_species, tempPlots, all.x=TRUE)
sampleDat = dcast(treeDat, plot_id + year_measured + latitude + longitude ~ id_spe, fill = 0, 
		value.var = "basal_area", fun.aggregate = function(x) as.integer(sum(x) > 0))
stateData = merge(sampleDat, tp_clim, all = 'T', by=c("plot_id", "year_measured"))


spName = '18032-ABI-BAL'

# to do:
# 1. read in the original random forest model (the one used to fit the STM
# 2. Project RF to the temporary plot grid
# 3. Project STM to temporary plot grid
# 4. for grid cells with multiple plots, select one at random
# 5. compute TSS for RF and for STM for both datasets

# 6. Read in original calibration P/A dataset
# 7. Aggregate PA to 1 or 0 for each climate grid cell for calibration and temp plot data
# 8. Build a new random forest on these aggregate data for the calibration set
# 9. Project RF and STM to the aggregated temporary plot dataset
# 10. Compute TSS for both models for the aggregated temp plot dataset