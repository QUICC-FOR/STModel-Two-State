#!/usr/bin/Rscript
# it's no longer necessary to run this interactively, since the predictor data is the same for all species
# this is not a model selection step, it is just selecting which variables aren't correlated with each other
# the SDM fitting is the model selection step, which is semi-automated (using step())
# note that since all species are the same, it's not really necessary to re-run this for each species
# should probably change this
if(!interactive()) stop(paste("SDM variable selection should be done interactively for each species",
	"This message means the script has likely changed and is being run via the makefile"
	"Run the '2_interactive_select_sdm_vars.r' in an interactive R session", sep="\n")

#-------------------
#  USER DEFINED VARIABLES
#  these must be set for each species
#-------------------
# spName = ""
# spName = "28731-ACE-SAC" # the script as written below is for this species
# spName = "18032-ABI-BAL"
spName = "28728-ACE-RUB"




# set seed - drawn from sample(1:1e6, 1)
set.seed(588533)
infile = paste("dat/", spName, "_processed.rdata", sep="")
load(infile)


# --------------------
#  variable selection 
# --------------------

# select vars to include in the SDM
# mean annual temp and total precip are included by default because we are using them as  
# the environmental variables in the CE model

# PCA
library(ade4)
var.pca = dudi.pca(stateData.subset[,climVarNames], scannf=FALSE, nf = 5)
var.pca$eig / sum(var.pca$eig)
varCor = cor(stateData.subset[,climVarNames])

contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

# procedure:
# examine nonCorVars at each step
# choose one with high inertia (in contrib) for the next step
# repeat until all non-correlated variables are included in selectedVars (below)
nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),
		"annual_mean_temp"]), names(varCor[which(abs(varCor[,"tot_annual_pp"])<0.7),
		"tot_annual_pp"]))
contrib[nonCorVars,]

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_seasonality"])<0.7),"pp_seasonality"]))
contrib[nonCorVars,]
nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_driest_period"])<0.7),"pp_driest_period"]))
contrib[nonCorVars,]
nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_warmest_quarter"])<0.7),"pp_warmest_quarter"]))
contrib[nonCorVars,]
nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"mean_diurnal_range"])<0.7),"mean_diurnal_range"]))
contrib[nonCorVars,]
nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"gdd_above_base_temp_period2"])<0.7),"gdd_above_base_temp_period2"]))
contrib[nonCorVars,]
selectedVars = c('annual_mean_temp', 'tot_annual_pp', 'pp_seasonality', 
		'mean_diurnal_range', 'pp_warmest_quarter', 'pp_driest_period', 
		'gdd_above_base_temp_period2', 'mean_temp_wettest_quarter')
		
sdmData = stateData.subset[,c('presence', selectedVars)]
save(sdmData, selectedVars, file=paste("dat/", spName, "_sdmData.rdata", sep=""))
