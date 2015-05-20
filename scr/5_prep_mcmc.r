#!/usr/bin/Rscript

library(argparse)
# handle command line arguments
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-f", "--fraction", default=0.25, type="double", help="fraction of data to use")
parser$add_argument("-r", "--rf", action="store_true", default=FALSE, help="random forest for prevalence (default: GLM)")
parser$add_argument("-m", "--gam", action="store_true", default=FALSE, help="use GAM for prevalence (default: GLM)")
argList = parser$parse_args()
spName = argList$species

datFile = paste("dat/", spName, "/", spName, "_transitions_projected.rds", sep="")
annealFile = paste("results/", spName, "/", spName, "_anneal_parameters.rds", sep="")
dat = readRDS(datFile)
annealPars = readRDS(annealFile)

cubic = (length(annealPars$par) > 10)
if(cubic)
{
	parNames = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
} else
{
	parNames = c('g0','g1','g2','g3','g4','e0','e1','e2','e3','e4')
}
mcmcDat = with(dat, data.frame(
	initial = state1,
	final = state2,
	env1 = annual_mean_temp,
	env2 = tot_annual_pp,
	interval = year2-year1))
	
mcmcDat = within(mcmcDat, 
{
	if(argList$gam)
	{
		prevalence1 = dat$expectedGAM
	} else if(argList$rf) {
		prevalence1 = dat$expectedRF
	} else {
		prevalence1 = dat$expectedGLM
	}
})


sampVar = rep(0.1, length(parNames))

parInfo = data.frame(
	name = parNames,
	initialValue = annealPars$par,
	samplerVariance = sampVar,
	priorMean = 0,
	priorSD = 10000,
	priorDist = "Normal")

mcmcFile = paste("run2/", spName, "/trans.txt", sep="")
parFile = paste("run2/", spName, "/inits.txt", sep="")

if(argList$fraction < 1)
{
	ind = sample(nrow(mcmcDat), argList$fraction*nrow(mcmcDat))
	mcmcDat = mcmcDat[ind,]
}
	
write.csv(mcmcDat, mcmcFile, row.names=FALSE)
write.csv(parInfo, parFile, row.names=FALSE)