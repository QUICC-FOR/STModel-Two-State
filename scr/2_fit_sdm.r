#!/usr/bin/Rscript

library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-p", "--power", default=2, type="integer", help="max exponent for model fitting")
argList = parser$parse_args()
spName = argList$species

# set seed - drawn from sample(1:1e6, 1)
## set.seed(588533)
infile = paste("dat/", spName, "/", spName, "_processed.rdata", sep="")
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


# --------------------
#  fit models
# --------------------

make_glm_model = function(vars, pow = 2)
{
	as.formula(paste("presence ~ ", paste(sapply(vars, function(v) paste(v, ifelse(pow > 1, paste(c("", sapply(2:pow, 
			function(p) paste("I(", v, "^", p, ")", sep=""))), collapse = ' + '), ""))), collapse = " + "), sep=""))
}

make_gam_scope = function(vars, pow)
{ lapply(vars, function(v) as.formula(paste("~ 1 + ", v, " + s(", v, ", ", pow, ")", sep=""))) }


## fit models
library(gam)
library(randomForest)

glm.mod = step(glm(make_glm_model(selectedVars, argList$power), family=binomial, data = 
		sdmData), k = log(nrow(sdmData)))

gam.mod = step.gam(gam(presence ~ ., family=binomial, data = sdmData),
		scope = make_gam_scope(selectedVars, argList$power))

rf.mod = randomForest(as.factor(presence) ~ . , data = sdmData, ntree = 500)

save(glm.mod, gam.mod, rf.mod, selectedVars, sdmData, file=paste("results/", spName, "/", spName, "_sdm_models.rdata", sep=""))





