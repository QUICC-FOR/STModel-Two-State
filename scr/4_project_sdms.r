#!/usr/bin/Rscript

library(argparse)
library(gam)
library(randomForest)
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
argList = parser$parse_args()


# project the SDMs to the climate grid
load(paste("results/", argList$species, "_sdm_models.rdata", sep=""))
load(paste("dat/", argList$species, "_sdmData.rdata", sep=""))
climGrid = readRDS(paste("dat/", argList$species, "_climData_scaled.rds", sep=""))
climGrid$glm.predict = predict(glm.mod, newdata=climGrid[,selectedVars], type='response')
climGrid$gam.predict = predict(gam.mod, newdata=climGrid[,selectedVars], type='response')
climGrid$rf.predict = predict(rf.mod, newdata=climGrid[,selectedVars], type='prob')[,2]
saveRDS(climGrid, paste("dat/", argList$species, "_climGrid_projected.rds", sep=""))


# project to the transition data for the model fitting
load(paste("dat/", argList$species, "_processed.rdata", sep=""))
transitionData.projected = transitionData.scaled
transitionData.projected$expectedGLM = predict(glm.mod, newdata = transitionData.scaled, type='response')
transitionData.projected$expectedGAM = predict(gam.mod, newdata = transitionData.scaled, type='response')
transitionData.projected$expectedRF = predict(rf.mod, newdata = transitionData.scaled, type='prob')[,2]
saveRDS(transitionData.projected, paste("dat/", argList$species, "_transitions_projected.rds", sep=""))


