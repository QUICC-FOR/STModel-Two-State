#!/usr/bin/Rscript
library(argparse)
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-p", "--power", default=2, type="integer", help="max exponent for model fitting")
argList = parser$parse_args()



# dear future me:
# here are three functions for you to try to figure out, because functional programming
# is AWESOME
# bon chance!

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
load(paste("dat/", argList$species, "_sdmData.rdata", sep=""))

glm.mod = step(glm(make_glm_model(selectedVars, argList$power), family=binomial, data = 
		sdmData), k = log(nrow(sdmData)))

gam.mod = step.gam(gam(presence ~ ., family=binomial, data = sdmData),
		scope = make_gam_scope(selectedVars, argList$power))

rf.mod = randomForest(as.factor(presence) ~ . , data = sdmData, ntree = 500)

save(glm.mod, gam.mod, rf.mod, file=paste("results/", argList$species, "_sdm_models.rdata", sep=""))




